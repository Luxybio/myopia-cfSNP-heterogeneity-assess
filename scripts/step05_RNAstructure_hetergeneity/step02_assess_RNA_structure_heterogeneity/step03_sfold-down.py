#!/usr/bin/env python3
# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
"""Batch-download Sfold prediction results (zip + recursive extraction)."""

import argparse
import os
import re
import sys
import time
import zipfile
from urllib.parse import urljoin, urlparse

try:
    import requests
except ImportError:
    sys.exit("Missing dependency, please run: pip install requests")

# Default configuration (modify as needed)
BASE_DIR = r"<LOCAL_DATA_ROOT>\Sfold-seq-build"
JOBS_FILE = {
    "WT": os.path.join(BASE_DIR, "sfold_submit", "WT", "WT_jobs.txt"),
    "MT": os.path.join(BASE_DIR, "sfold_submit", "MT", "MT_jobs.txt"),
}
OUT_ROOT = os.path.join(BASE_DIR, "sfold-down")

TIMEOUT     = 120     # Timeout extended for potentially large zips
RETRY       = 3       # Retry count when zip is not ready
WAIT        = 6       # Seconds to wait between retries (server needs time to package)
DELAY       = 0.5     # Polite delay between tasks (seconds)
HEADERS     = {"User-Agent": "Mozilla/5.0 (compatible; sfold-down/1.0)"}

# Candidate files for per-file mirroring fallback (when zip fails; from the observed structure of old directory 15920)
FIXED_FILES = [
    "srna.html", "jobinfo.html", "bp.out",
    "output/fe.out", "output/pdf.out", "output/sclass.out", "output/cdf.out",
    "output/sample.out", "output/smfe.out", "output/sirna.out",
    "output/sirna_de.out", "output/sirna_f.out", "output/sirna_s.out",
    "output/oligo.out", "output/oligo_de.out", "output/oligo_f.out",
    "output/Dharmacon_thermo.out",
    "cdf/cdf.out", "cdf/cdf.pdf",
    "circles/ccentroids.pdf", "circles/ec_mfe.pdf",
    "stabprof/stability.out",
    "probprof/complete.pdf", "probprof/loopr.out", "probprof/sstrand.out",
    "probprof/bplot_1.pdf", "probprof/extplot_1.pdf", "probprof/hplot_1.pdf",
    "probprof/iplot_1.pdf", "probprof/mplot_1.pdf", "probprof/regional_1.pdf",
    "structures/ecentroid.ct", "structures/ecentroid.pdf",
    "structures/mfold.ct", "structures/mfold.pdf",
    "mds/mds.out", "triplet/triplet.out",
    "2dhist/2dhist.out",
]
PATTERN_FILES = [
    ("2dhist/c{i:02d}.2dhist.out", 20),
    ("structures/c{i:02d}.ccentroid.ct", 20),
    ("structures/c{i:02d}.ccentroid.pdf", 20),
    ("structures/structure{i}.out", 20),
    ("structures/structure{i}.pdf", 20),
]
MISS_TOL = 2     # Stop a pattern family after this many consecutive 404s


# Basic utilities
def norm_jobs_arg(path):
    """Allow passing MT_jobs / WT_jobs without extension, automatically appending .txt etc."""
    if os.path.isfile(path):
        return path
    for ext in (".txt", ".tsv", ".list"):
        if os.path.isfile(path + ext):
            return path + ext
    return path


def read_jobs(path):
    """Read jobs file -> [(jobname, url, header), ...]"""
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            jobname, url = parts[0].strip(), parts[1].strip()
            header = parts[2].strip() if len(parts) > 2 else ""
            if jobname and url.startswith("http"):
                rows.append((jobname, url, header))
    return rows


def job_base(url):
    """.../output/<dt>.<pid>/srna.html -> (base_url, dt, pid)"""
    p = urlparse(url)
    m = re.match(r"^(.*/output/(\d+)\.(\d+))/", p.path)
    if not m:
        return None, None, None
    return "%s://%s%s/" % (p.scheme, p.netloc, m.group(1)), m.group(2), m.group(3)


def save_bytes(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".part"
    with open(tmp, "wb") as fh:
        fh.write(content)
    os.replace(tmp, path)


def is_zip_bytes(b):
    return b is not None and len(b) > 4 and b[:2] == b"PK"


_EXTRACTED = set()   # Absolute paths of already extracted zips (cross-call deduplication)


def extract_zip_recursive(zip_path, dest_dir, log=print):
    """
    Extract zip; continue extracting any nested zips found in the extracted contents (recursive).
    Intermediate zips are deleted after extraction; the outermost user-downloaded zip is kept.
    Returns the total number of extracted entries.
    """
    n_entries = 0
    while True:
        todo = None
        for root, _d, files in os.walk(dest_dir):
            for f in files:
                if f.lower().endswith(".zip"):
                    p = os.path.join(root, f)
                    if p not in _EXTRACTED:
                        todo = p
                        break
            if todo:
                break
        if todo is None:
            break
        _EXTRACTED.add(todo)
        try:
            with zipfile.ZipFile(todo) as zf:
                names = zf.namelist()
                zf.extractall(dest_dir)
            n_entries += len(names)
            log("    Extracted %s (%d entries)" % (os.path.relpath(todo, dest_dir), len(names)))
        except zipfile.BadZipFile:
            log("    [WARN] Invalid zip: %s" % todo)
            continue
        if todo != zip_path:      # Delete intermediate zips after use
            try:
                os.remove(todo)
            except OSError:
                pass
    return n_entries


# Zip download method (corresponds to the two Download buttons on the page)
def download_general_zip(sess, base, dt, pid, jdir, log=print):
    """General output: POST download.pl -> follow two meta redirects to download <pid>.zip."""
    try:
        r = sess.post("https://sfold.wadsworth.org/cgi-bin/download.pl",
                      data={"dt": dt, "pid": pid, "triplet": "GUC", "ftype": "0"},
                      timeout=TIMEOUT)
    except requests.RequestException:
        return None
    if r.status_code != 200:
        return None

    # Redirect chain: download.pl -> downloading.pl (triggers packaging) -> /output/<dt>.<pid>/<pid>.zip
    zurl = urljoin("https://sfold.wadsworth.org/", "cgi-bin/downloading.pl?dt=%s&pid=%s&ftype=0"
                   % (dt, pid))
    m = re.search(r'url=([^">]+)', r.text)
    if m:
        zurl = urljoin("https://sfold.wadsworth.org/", m.group(1).strip())

    for attempt in range(1, RETRY + 1):
        try:
            r2 = sess.get(zurl, timeout=TIMEOUT)
        except requests.RequestException:
            r2 = None
        if r2 is not None and r2.status_code == 200 and is_zip_bytes(r2.content):
            zpath = os.path.join(jdir, "%s.zip" % pid)
            save_bytes(zpath, r2.content)
            return zpath
        # HTML page: keep parsing meta refresh (downloading.pl -> <pid>.zip)
        if r2 is not None and r2.status_code == 200:
            m2 = re.search(r'url=([^">]+)', r2.text)
            if m2:
                zurl = urljoin(zurl, m2.group(1).strip())
        log("    zip not ready (%s), retrying in %ds %d/%d"
            % (getattr(r2, "status_code", "network error"), WAIT, attempt, RETRY))
        time.sleep(WAIT)
    return None


# Fallback: per-file mirroring
def mirror_files(sess, base, jdir, force=False):
    ok = miss = 0

    def fetch(rel):
        nonlocal ok, miss
        dst = os.path.join(jdir, rel)
        if os.path.exists(dst) and os.path.getsize(dst) > 0 and not force:
            ok += 1
            return True
        try:
            r = sess.get(urljoin(base, rel), timeout=TIMEOUT)
        except requests.RequestException:
            miss += 1
            return False
        if r.status_code == 200 and r.content:
            save_bytes(dst, r.content)
            ok += 1
            return True
        miss += 1
        return False

    for rel in FIXED_FILES:
        fetch(rel)
    for tpl, limit in PATTERN_FILES:
        streak = 0
        for i in range(1, limit + 1):
            hit = fetch(tpl.format(i=i))
            streak = 0 if hit else streak + 1
            if streak >= MISS_TOL:
                break
    return ok, miss


# Single task
def download_job(sess, jobname, url, out_dir, force=False, log=print):
    base, dt, pid = job_base(url)
    if base is None:
        return {"job": jobname, "mode": "BAD_URL", "files": 0, "zips": ""}

    jdir = os.path.join(out_dir, jobname)
    os.makedirs(jdir, exist_ok=True)
    # Completion marker (resume support): stabprof/stability.out exists only inside the general zip,
    # and the MIRROR fallback cannot retrieve it, so tasks with only mirroring are treated as incomplete.
    marker = os.path.join(jdir, pid, "stabprof", "stability.out")
    if os.path.exists(marker) and not force:
        return {"job": jobname, "mode": "SKIP", "files": 0, "zips": ""}

    zips = []
    n = 0
    zp = download_general_zip(sess, base, dt, pid, jdir, log=log)
    if zp:
        n += extract_zip_recursive(zp, jdir, log=log)
        zips.append("%s.zip" % pid)

    mode = "ZIP"
    if not zp:                                   # Zip failed -> per-file fallback
        mode = "MIRROR"
        ok, _miss = mirror_files(sess, base, os.path.join(jdir, pid), force=force)
        n += ok

    return {"job": jobname, "mode": mode, "files": n, "zips": "+".join(zips)}


# Main workflow
def run(jobs_path, out_root, force=False, only=None):
    jobs_path = norm_jobs_arg(jobs_path)
    if not os.path.isfile(jobs_path):
        sys.exit("[FATAL] jobs file does not exist: %s" % jobs_path)
    stype = os.path.basename(jobs_path).split("_")[0] or "JOB"
    out_dir = os.path.join(out_root, only or stype)
    os.makedirs(out_dir, exist_ok=True)

    rows = read_jobs(jobs_path)
    if not rows:
        sys.exit("[FATAL] no usable rows in jobs file: %s" % jobs_path)
    print("[INFO] %s: %d tasks -> %s" % (os.path.basename(jobs_path), len(rows), out_dir))

    sess = requests.Session()
    sess.headers.update(HEADERS)

    stat = []
    for i, (jobname, url, _h) in enumerate(rows, 1):
        print("[%d/%d] %s" % (i, len(rows), jobname), flush=True)
        res = download_job(sess, jobname, url, out_dir, force=force)
        print("    mode=%s entries=%d zip=%s" % (res["mode"], res["files"], res["zips"] or "-"),
              flush=True)
        stat.append(res)
        time.sleep(DELAY)

    sumpath = os.path.join(out_dir, "download_summary.txt")
    with open(sumpath, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("jobname\tmode\tentries\tzips\n")
        for r in stat:
            fh.write("%s\t%s\t%d\t%s\n" % (r["job"], r["mode"], r["files"], r["zips"]))
    print("\n===== Summary =====")
    print("Tasks %d; ZIP %d, MIRROR fallback %d, skipped %d"
          % (len(stat),
             sum(r["mode"] == "ZIP" for r in stat),
             sum(r["mode"] == "MIRROR" for r in stat),
             sum(r["mode"] == "SKIP" for r in stat)))
    print("Summary table: %s" % sumpath)


def main():
    ap = argparse.ArgumentParser(description="Batch download Sfold results (zip + recursive extraction)")
    ap.add_argument("--jobs", default="", help="Path to jobs list (default uses --type to get the default path)")
    ap.add_argument("--type", choices=["WT", "MT", "BOTH"], default="BOTH")
    ap.add_argument("--out", default=OUT_ROOT, help="Output root directory")
    ap.add_argument("--force", action="store_true", help="Ignore completion marker and re-download")
    args = ap.parse_args()

    if args.jobs:
        stype = os.path.basename(norm_jobs_arg(args.jobs)).split("_")[0]
        run(args.jobs, args.out, force=args.force,
            only=stype if stype in ("WT", "MT") else None)
    else:
        types = ["WT", "MT"] if args.type == "BOTH" else [args.type]
        for t in types:
            run(JOBS_FILE[t], args.out, force=args.force, only=t)


if __name__ == "__main__":
    main()
