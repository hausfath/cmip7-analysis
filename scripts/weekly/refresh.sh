#!/bin/bash
# Weekly CMIP7 refresh, run by launchd from a dedicated clone of this repository.
#
#   scan ESGF -> nothing new: stop (no commit)
#             -> new files: acquire, rerun steps 4-8, mechanical checks (check.py),
#                Claude updates README/METHODS and gives a verdict,
#                clean  -> commit and push to main (GitHub Pages updates)
#                review -> commit to branch weekly-refresh and open/update a pull request
#
# Layout: $CMIP7_AUTO_DIR (default ~/cmip7-weekly)/repo = clone with its own .venv,
#         $CMIP7_AUTO_DIR/logs/<date>/ = run.log, scan.json, check.json, verdict.json, summary.md
# Flags: --dry-run  stop after the mechanical checks (no Claude step, no commit, no push)
set -uo pipefail

main() {
  local AUTO="${CMIP7_AUTO_DIR:-$HOME/cmip7-weekly}"
  local REPO="$AUTO/repo" BR=weekly-refresh DRY=0
  [ "${1:-}" = "--dry-run" ] && DRY=1
  local PY="$REPO/.venv/bin/python"
  local DATE; DATE=$(date +%F)
  local LOG="$AUTO/logs/$DATE"
  mkdir -p "$LOG"
  exec >>"$LOG/run.log" 2>&1
  export PATH="/opt/homebrew/bin:$HOME/.local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
  echo "=== $(date '+%F %T') start (dry-run=$DRY)"

  if ! mkdir "$AUTO/lock" 2>/dev/null; then echo "another run holds $AUTO/lock; exiting"; return 0; fi
  trap 'rmdir "$AUTO/lock"' EXIT
  caffeinate -i -w $$ &      # keep the Mac awake while this runs

  notify() { osascript -e "display notification \"$1\" with title \"CMIP7 weekly refresh\"" >/dev/null 2>&1 || true
             printf '%s\t%s\n' "$DATE" "$1" >>"$AUTO/logs/history.tsv"; }
  fail() { echo "FAILED: $1"; notify "Failed: $1 (log: $LOG/run.log)"; exit 1; }

  cd "$REPO" || fail "no clone at $REPO"
  git fetch -q origin || fail "git fetch"
  # An open refresh PR is built on, so the next run does not redo or bypass unreviewed work
  local OPEN_PR; OPEN_PR=$(gh pr list --head "$BR" --state open --json number -q '.[0].number' 2>/dev/null)
  git reset -q --hard && git clean -fdq
  if [ -n "$OPEN_PR" ]; then
    git checkout -q -B "$BR" "origin/$BR" && git merge -q --no-edit origin/main || fail "merge main into $BR"
  else
    git checkout -q -B "$BR" origin/main || fail "checkout"
  fi

  "$PY" -I scripts/build_manifest.py || fail "manifest query (ESGF index unreachable?)"
  "$PY" -I scripts/weekly/scan.py "$LOG"; local rc=$?
  if [ $rc -eq 0 ]; then
    git reset -q --hard && git clean -fdq
    echo "nothing new on ESGF"; printf '%s\t%s\n' "$DATE" "nothing new" >>"$AUTO/logs/history.tsv"
    return 0
  fi
  [ $rc -eq 10 ] || fail "scan.py"

  while read -r f; do [ -n "$f" ] && git rm -q --ignore-unmatch "$f"; done <"$LOG/refetch.txt"
  # Documented exception (METHODS §1.1): TLS checks off for the CAS node only while its
  # certificate is invalid (curl exit 60); acquire still requires an exact sha256 match there.
  local INSECURE=()
  curl -s -o /dev/null --max-time 30 https://esg-earthlab.iap.ac.cn/ ; [ $? -eq 60 ] && INSECURE=(--insecure-host esg-earthlab.iap.ac.cn)
  local STEPFAIL=""
  for pass in 1 2; do   # acquisition is resumable; a second pass retries transient failures
    "$PY" -I scripts/acquire_global_means.py ${INSECURE[@]+"${INSECURE[@]}"} > "$LOG/acquire_$pass.txt" 2>&1
    tail -3 "$LOG/acquire_$pass.txt"
    "$PY" -I scripts/status.py | grep -q "all complete" && break
  done
  for s in "sensitivity.py cmip7" gmst.py figures_gmst.py figures_sensitivity.py dataset_table.py build_report.py; do
    echo "--- $s"; "$PY" -I scripts/$s > "$LOG/step_${s%% *}.txt" 2>&1 || { STEPFAIL="scripts/$s failed (see $LOG/step_${s%% *}.txt)"; break; }
  done
  "$PY" -I scripts/weekly/check.py "$LOG" || STEPFAIL="${STEPFAIL:-check.py failed}"
  if [ $DRY -eq 1 ]; then echo "dry run: stopping before the Claude step"; return 0; fi

  local VERDICT=review
  if [ -z "$STEPFAIL" ]; then
    git diff -- outputs/ > "$LOG/outputs.diff"; git diff --stat > "$LOG/diffstat.txt"
    sed -e "s|{LOG}|$LOG|g" -e "s|{DATE}|$DATE|g" scripts/weekly/claude_prompt.md > "$LOG/prompt.md"
    claude -p "$(cat "$LOG/prompt.md")" --model opus --permission-mode dontAsk \
      --allowedTools "Read,Glob,Grep,Edit(./README.md),Edit(./METHODS.md),Edit(/$LOG/verdict.json),Edit(/$LOG/summary.md),Bash(git diff *),Bash(git show *),Bash(git log *),Bash(git status)" \
      > "$LOG/claude.txt" 2>&1 || echo "claude exited non-zero"
    VERDICT=$("$PY" -I -c "import json,sys; print(json.load(open(sys.argv[1]))['verdict'])" "$LOG/verdict.json" 2>/dev/null || echo review)
    "$PY" -I -c "import json,sys; sys.exit(0 if json.load(open(sys.argv[1]))['clean'] else 1)" "$LOG/check.json" || VERDICT=review
  fi
  [ -n "$OPEN_PR" ] && VERDICT=review   # never publish on top of an unreviewed refresh
  [ -s "$LOG/summary.md" ] || { echo "Weekly refresh $DATE."; [ -n "$STEPFAIL" ] && echo "$STEPFAIL"; } > "$LOG/summary.md"
  echo "verdict: $VERDICT"

  git add -A
  { echo "Weekly refresh: ESGF snapshot $DATE ($VERDICT)"; echo; cat "$LOG/summary.md"
    [ -n "$STEPFAIL" ] && { echo; echo "Pipeline failure: $STEPFAIL"; }
    echo; echo "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"; } > "$LOG/commit_msg.txt"
  git commit -q -F "$LOG/commit_msg.txt" || fail "git commit"

  if [ "$VERDICT" = clean ] && git push -q origin HEAD:main; then
    notify "Published: ESGF snapshot $DATE pushed to main"
    return 0
  fi
  if [ -n "$OPEN_PR" ]; then git push -q origin "HEAD:refs/heads/$BR"
  else git push -q --force origin "HEAD:refs/heads/$BR"; fi || fail "push $BR"   # no open PR: any old $BR is merged or abandoned
  { cat "$LOG/summary.md"; echo; echo "Logs on the automation machine: \`$LOG\`"
    echo; echo "🤖 Generated with [Claude Code](https://claude.com/claude-code)"; } > "$LOG/pr_body.md"
  if [ -n "$OPEN_PR" ]; then
    gh pr comment "$OPEN_PR" --body-file "$LOG/pr_body.md" >/dev/null
  else
    gh pr create --base main --head "$BR" --title "Weekly CMIP7 refresh needs review ($DATE)" --body-file "$LOG/pr_body.md" >/dev/null
  fi
  notify "Needs review: pull request on $BR ($DATE)"
}
main "$@"; exit $?
