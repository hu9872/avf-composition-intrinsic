"""多线程分块下载器：应对 NCBI 单连接限速（实测 ~40 KB/s，16 线程可聚合）。"""
import urllib.request, os, sys, time, hashlib, socket
from concurrent.futures import ThreadPoolExecutor, as_completed

socket.setdefaulttimeout(180)

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/{grp}/{acc}/suppl/{f}"
OUT = r"D:/AVF_GSE_2611/data/raw"
TMP = os.path.join(OUT, "_chunks")
CHUNK = 4 << 20
THREADS = 12

TASKS = [
    ("GSE305423", "GSE305423_RAW.tar"),
    ("GSE220796", "GSE220796_rVazquez-Padron_counts.txt.gz"),
    ("GSE119296", "GSE119296_Vazquez-Padron_RNAseq_080718_counts.txt.gz"),
    ("GSE314264", "GSE314264_Final_counts-AVF_timepoints.csv.gz"),
    ("GSE39488", "GSE39488_RAW.tar"),
]


def url(acc, f):
    return BASE.format(grp=acc[:-3] + "nnn", acc=acc, f=f)


def head(acc, f):
    for _ in range(3):
        try:
            r = urllib.request.urlopen(urllib.request.Request(url(acc, f), method="HEAD",
                                                              headers={"User-Agent": "Mozilla/5.0"}))
            return int(r.headers["Content-Length"])
        except Exception as e:
            print("  head retry:", type(e).__name__, flush=True)
            time.sleep(3)
    return None


def fetch(acc, f, i, total, dst):
    if os.path.exists(dst) and os.path.getsize(dst) == min(CHUNK, total - i * CHUNK):
        return dst
    s = i * CHUNK
    e = min(s + CHUNK, total) - 1
    for att in range(6):
        try:
            time.sleep(0.15 * (i % 7))
            h = {"User-Agent": "Mozilla/5.0", "Range": f"bytes={s}-{e}"}
            d = urllib.request.urlopen(urllib.request.Request(url(acc, f), headers=h)).read()
            if len(d) == e - s + 1:
                with open(dst, "wb") as fh:
                    fh.write(d)
                return dst
            raise IOError(f"short read {len(d)}/{e-s+1}")
        except Exception as ex:
            code = getattr(ex, "code", None)
            if att == 5:
                raise
            wait = 5 * (2 ** att) if code in (429, 503, 500) else 2 * (att + 1)
            time.sleep(min(wait, 90))


def get(acc, f):
    final = os.path.join(OUT, f)
    if os.path.exists(final) and os.path.getsize(final) > 0:
        print(f"[skip] {f}", flush=True)
        return True
    total = head(acc, f)
    if not total:
        print(f"[FAIL-head] {f}", flush=True)
        return False
    nd = (total + CHUNK - 1) // CHUNK
    d = os.path.join(TMP, f"{acc}_{f}")
    os.makedirs(d, exist_ok=True)
    print(f"== {f}  {total/1e6:.1f} MB  {nd} chunks", flush=True)
    t0 = time.time()
    done = 0
    for rnd in (1, 2):
        missing = [i for i in range(nd) if not os.path.exists(os.path.join(d, f"{i:05d}.part"))]
        if not missing:
            break
        print(f"   pass{rnd}: {len(missing)} chunks to fetch", flush=True)
        with ThreadPoolExecutor(max_workers=THREADS) as ex:
            futs = {ex.submit(fetch, acc, f, i, total, os.path.join(d, f"{i:05d}.part")): i
                    for i in missing}
            for fu in as_completed(futs):
                try:
                    fu.result()
                except Exception as ex_:
                    print(f"   chunk {futs[fu]} failed: {type(ex_).__name__}", flush=True)
                done += 1
                if done % 25 == 0:
                    got = sum(os.path.getsize(os.path.join(d, x)) for x in os.listdir(d))
                    print(f"   {got/1e6:.0f}/{total/1e6:.0f} MB  "
                          f"{got/1e3/max(time.time()-t0,1):.0f} KB/s", flush=True)
        time.sleep(3)
    got_files = [i for i in range(nd) if os.path.exists(os.path.join(d, f"{i:05d}.part"))]
    if len(got_files) < nd:
        print(f"[INCOMPLETE] {f}: {len(got_files)}/{nd} chunks", flush=True)
        return False
    with open(final, "wb") as out:
        for i in range(nd):
            with open(os.path.join(d, f"{i:05d}.part"), "rb") as fh:
                while True:
                    b = fh.read(1 << 22)
                    if not b:
                        break
                    out.write(b)
    sz = os.path.getsize(final)
    ok = sz == total
    print(f"[{'ok' if ok else 'BAD'}] {f} {sz/1e6:.1f} MB", flush=True)
    if ok:
        h = hashlib.md5()
        with open(final, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 22), b""):
                h.update(b)
        with open(os.path.join(OUT, "md5.txt"), "a") as fh:
            fh.write(f"{h.hexdigest()}  {f}\n")
        for x in os.listdir(d):
            os.remove(os.path.join(d, x))
        os.rmdir(d)
    return ok


if __name__ == "__main__":
    only = sys.argv[1:] or None
    for acc, f in TASKS:
        if only and f not in only and acc not in only:
            continue
        print("==", acc, f, flush=True)
        get(acc, f)
    print("ALL DONE", flush=True)
