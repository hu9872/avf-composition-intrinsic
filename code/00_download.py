import urllib.request, os, sys, time, socket

socket.setdefaulttimeout(120)
OUT = r"D:/AVF_GSE_2611/data/raw"
os.makedirs(OUT, exist_ok=True)

TASKS = [
    ("GSE305423", "GSE305423_RAW.tar"),
    ("GSE220796", "GSE220796_rVazquez-Padron_counts.txt.gz"),
    ("GSE119296", "GSE119296_Vazquez-Padron_RNAseq_080718_counts.txt.gz"),
    ("GSE314264", "GSE314264_Final_counts-AVF_timepoints.csv.gz"),
    ("GSE244317", None),
    ("GSE39488", "GSE39488_RAW.tar"),
]

def url_for(acc, f):
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{f}"

def get(acc, f):
    dst = os.path.join(OUT, f)
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        print(f"[skip] {f} exists {os.path.getsize(dst)/1e6:.1f} MB", flush=True)
        return True
    tmp = dst + ".part"
    for attempt in (1, 2):
        try:
            req = urllib.request.Request(url_for(acc, f), headers={"User-Agent": "Mozilla/5.0"})
            r = urllib.request.urlopen(req, timeout=120)
            total = int(r.headers.get("Content-Length", 0))
            done = 0
            t0 = time.time()
            with open(tmp, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
                    done += len(chunk)
                    if total and done % (50 << 20) < (1 << 20):
                        pct = 100 * done / total
                        sp = done / 1e6 / max(time.time() - t0, 1)
                        print(f"  {f}: {pct:5.1f}%  {done/1e6:.0f}/{total/1e6:.0f} MB  {sp:.1f} MB/s", flush=True)
            os.replace(tmp, dst)
            print(f"[ok] {f} {os.path.getsize(dst)/1e6:.1f} MB", flush=True)
            return True
        except Exception as e:
            print(f"[fail attempt {attempt}] {f}: {type(e).__name__} {str(e)[:120]}", flush=True)
            time.sleep(5)
    return False

if __name__ == "__main__":
    for acc, f in TASKS:
        if f is None:
            continue
        print(f"== {acc} / {f}", flush=True)
        get(acc, f)
    print("ALL DONE", flush=True)
