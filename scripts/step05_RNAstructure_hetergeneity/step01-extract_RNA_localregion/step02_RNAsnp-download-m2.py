# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
"""Batch-download RNAsnp Mode 2 result files from recorded URLs (requests only)."""

import os
import re
import traceback

import requests

element = "LncRNA"   # Element name (must match submission script)
seqlen = "over"

alldir = "<LOCAL_DATA_ROOT>/negative-RNAsnp-rerun/"
result_file = alldir + seqlen + "1wan_" + element + "_RNAsnp_result.txt"
out_base = alldir + element + "-" + seqlen + "-RNAsnp/"

# Download targets: subdirectory -> (link regex pattern, saved extension)
# Patterns are matched against hrefs on the result page with re.search; {job} is the jobid
DOWNLOAD_TARGETS = {
    "interval-p":       (r"RNAsnp_results_{job}\.csv$", ".txt"),
    "wt-mt-prob-pdf":   (r"_dpp\.pdf$", ".pdf"),
    "wt-mt-prob-value": (r"_dpp\.txt$", ".txt"),
    "wt-structure-pdf": (r"_WT_dp\.pdf$", ".pdf"),
    "wt-structure-dp":  (r"_WT_dp\.txt$", ".txt"),
    "mt-structure-pdf": (r"_MT_dp\.pdf$", ".pdf"),
    "mt-structure-dp":  (r"_MT_dp\.txt$", ".txt"),
}

# Expected content magic bytes for each file type (download validation to avoid saving HTML error pages as result files)
MAGIC = {".pdf": b"%PDF", ".txt": None}   # txt is not strictly validated (dp/value are numeric matrices)


def load_result_lines():
    """Read the submission result file and return a list of (sample_name, result_url)."""
    records = []
    with open(result_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip("\n")
            if not line:
                continue
            info = line.split("\t")
            if len(info) >= 3:
                records.append((info[1], info[2]))
    return records


def abs_link(href):
    if href.startswith("http"):
        return href
    return "https://rth.dk" + href


def download_all(sample, url, session):
    """Visit one result page and download all 7 file types."""
    resp = session.get(url, timeout=60)
    resp.raise_for_status()
    html = resp.text

    m = re.search(r"jobid=([0-9\-]+)", url)
    job = m.group(1) if m else ""

    # Collect all download links on the result page
    hrefs = re.findall(r'href="([^"]+)"', html)

    ok = True
    for sub, (pattern, ext) in DOWNLOAD_TARGETS.items():
        out_dir = out_base + sub + "/"
        os.makedirs(out_dir, exist_ok=True)
        out_path = out_dir + sample + ext

        if os.path.isfile(out_path):
            print(f"  [Skip] {out_path} already exists")
            continue

        pat = re.compile(pattern.format(job=re.escape(job)))
        link = None
        for h in hrefs:
            path = h.split("?")[0]
            if pat.search(path):
                link = abs_link(h)
                break
        if link is None:
            print(f"  [Warning] {sub}: no matching link found on result page ({pattern}), skipping")
            ok = False
            continue

        r2 = session.get(link, timeout=120)
        content = r2.content
        # Validate: pdf must start with %PDF; txt must not be HTML/ZIP
        if ext == ".pdf" and not content.startswith(b"%PDF"):
            print(f"  [Warning] {sub}: downloaded content is not a PDF, skipping save")
            ok = False
            continue
        if ext == ".txt" and (content[:2] == b"PK" or content[:1] == b"<"):
            print(f"  [Warning] {sub}: downloaded content is ZIP/HTML instead of a numeric table, skipping save")
            ok = False
            continue

        with open(out_path, "wb") as f:
            f.write(content)
        print(f"  [Downloaded] {out_path} ({len(content)} bytes)")
    return ok


def main():
    if not os.path.isfile(result_file):
        raise SystemExit(f"Result file does not exist: {result_file}")

    records = load_result_lines()
    print(f"{len(records)} records to process")

    session = requests.Session()
    failed = []

    for sample, url in records:
        # Skip failed submissions (result file may contain invalid URLs like submit.php#input)
        if "results.php?jobid=" not in url:
            print(f"[Skip] Invalid result page URL (submission may have failed): {sample}\t{url}")
            continue
        print(f"{sample}\t{url}")
        try:
            if not download_all(sample, url, session):
                failed.append(sample)
        except Exception:
            print(f"[Failed] {sample} download exception, see error below")
            traceback.print_exc()
            failed.append(sample)

    print(f"\nDone. Successful: {len(records) - len(failed)} | with missing items/failed: {len(failed)}")
    if failed:
        print("Failed sample list:")
        for s in failed:
            print("  -", s)


if __name__ == "__main__":
    main()
