import json, sys, time, threading, pandas as pd, requests
m = pd.read_csv("inventory/manifest_latest.csv")
def pull(url, secs, out):
    t0=time.time(); n=0
    with requests.get(url, stream=True, timeout=60) as r:
        for ch in r.iter_content(1<<20):
            n+=len(ch)
            if time.time()-t0>secs: break
    out.append(n/(time.time()-t0))
for src, npar in [("CanESM6-0-MR",1),("CanESM6-0-MR",4),("UKCM2a-0-HH",1),("UKCM2a-0-HH",4),("ACCESS-ESM1-6",4)]:
    urls=[json.loads(u)[0] for u in m[(m.source_id==src)&(m["size"]>8e7)].urls.head(npar)]
    out=[]; th=[threading.Thread(target=pull,args=(u,20,out)) for u in urls]
    [t.start() for t in th]; [t.join() for t in th]
    print(src, npar, "streams: total MB/s", round(sum(out)/1e6,2), [round(x/1e6,2) for x in out], flush=True)
