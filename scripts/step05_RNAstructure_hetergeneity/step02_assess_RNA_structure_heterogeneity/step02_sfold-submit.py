#!/usr/bin/env python3
# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
"""Batch-submit WT/MT local sequences to the Sfold online server."""

import argparse
import os
import time

from selenium import webdriver
from selenium.webdriver.common.by import By

# Paths and parameters (modify here as needed)
BASE_DIR   = r"<LOCAL_DATA_ROOT>\Sfold-seq-build"
TYPE_SEQ_DIR = {                     # type -> sequence folder
    "WT": os.path.join(BASE_DIR, "WT"),
    "MT": os.path.join(BASE_DIR, "MT"),
}
TYPE_OUT_DIR = {                     # type -> job list output folder
    "WT": os.path.join(BASE_DIR, "sfold_submit", "WT"),
    "MT": os.path.join(BASE_DIR, "sfold_submit", "MT"),
}

SFOLD_URL     = "https://sfold.wadsworth.org/cgi-bin/index.pl"
SUBMIT_WAIT   = 55      # Seconds to wait after submission for the server to generate results (old script used 55s)
AFTER_WAIT    = 10      # Seconds to wait after recording the URL
MAX_RETRY     = 2       # Retry count for a failed single task
CHROMEDRIVER_PATH = ""  # Leave empty to let Selenium Manager auto-manage (recommended for selenium>=4.6)

# Page elements (XPath verified with the old script)
XPATH_SRNA_ICON  = "/html/body/center/table/tbody/tr[3]/td/table[1]/tbody/tr[3]/td[2]/table/tbody/tr[4]/td[1]/a/img"
XPATH_JOB_INPUT  = "/html/body/center/table/tbody/tr[3]/td/form/table/tbody/tr[4]/td[2]/input"
XPATH_SEQ_INPUT  = "/html/body/center/table/tbody/tr[3]/td/form/table/tbody/tr[5]/td[2]/textarea"
XPATH_SUBMIT_BTN = "/html/body/center/table/tbody/tr[3]/td/form/table/tbody/tr[12]/td[2]/input[1]"


def read_seq(path):
    """Read a single FASTA entry, return (header, seq)."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        raw = fh.read()
    header, buf = None, []
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip()
        else:
            buf.append(line)
    return header or os.path.basename(path), "".join(buf)


def load_done(jobfile):
    """Read already submitted jobnames; supports resume."""
    done = set()
    if os.path.exists(jobfile):
        with open(jobfile, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.split("\t")
                if parts and parts[0].strip():
                    done.add(parts[0].strip())
    return done


def make_driver():
    if CHROMEDRIVER_PATH:
        from selenium.webdriver.chrome.service import Service
        return webdriver.Chrome(service=Service(CHROMEDRIVER_PATH))
    return webdriver.Chrome()   # Selenium >= 4.6 downloads/manages the driver automatically


def submit_one(driver, jobname, seqfa):
    """Submit one sequence and return the result page URL; raises exception on failure."""
    driver.get(SFOLD_URL)

    driver.find_element(By.XPATH, XPATH_SRNA_ICON).click()      # Enter sRNA folding entry
    driver.find_element(By.XPATH, XPATH_JOB_INPUT).send_keys(jobname)
    driver.find_element(By.XPATH, XPATH_SEQ_INPUT).send_keys(seqfa)
    driver.find_element(By.XPATH, XPATH_SUBMIT_BTN).click()

    time.sleep(SUBMIT_WAIT)          # Wait for the server to generate results
    return driver.current_url


def run_type(stype, dry_run=False):
    seqdir  = TYPE_SEQ_DIR[stype]
    outdir  = TYPE_OUT_DIR[stype]
    os.makedirs(outdir, exist_ok=True)
    jobfile = os.path.join(outdir, "%s_jobs.txt" % stype)

    files = sorted(f for f in os.listdir(seqdir) if f.lower().endswith(".txt"))
    done  = load_done(jobfile)
    todo  = [f for f in files if "%s-%s" % (stype, f[:-4]) not in done]

    print("[%s] %d sequence files, %d completed, %d remaining"
          % (stype, len(files), len(files) - len(todo), len(todo)))
    if dry_run:
        for f in todo:
            print("   would submit:", "%s-%s" % (stype, f[:-4]))
        return 0, 0

    ok = fail = 0
    driver = make_driver()
    try:
        for f in todo:
            jobname = "%s-%s" % (stype, f[:-4])
            header, seqfa = read_seq(os.path.join(seqdir, f))
            print("  submitting:", jobname, "(len=%d)" % len(seqfa), flush=True)

            success = False
            for attempt in range(1, MAX_RETRY + 1):
                try:
                    url = submit_one(driver, jobname, seqfa)
                    with open(jobfile, "a", encoding="utf-8", newline="\n") as out:
                        out.write("%s\t%s\t%s\n" % (jobname, url, header))
                    print("    OK ->", url, flush=True)
                    ok += 1
                    success = True
                    break
                except Exception as exc:
                    print("    attempt %d/%d failed: %s" % (attempt, MAX_RETRY, exc), flush=True)
                    time.sleep(5)
                    try:                     # Page state may be broken; rebuild browser
                        driver.quit()
                    except Exception:
                        pass
                    driver = make_driver()
            if not success:
                fail += 1
                print("    FAILED:", jobname, flush=True)

            time.sleep(AFTER_WAIT)
    finally:
        try:
            driver.quit()
        except Exception:
            pass

    print("[%s] Done: success %d, failed %d, job list: %s" % (stype, ok, fail, jobfile))
    return ok, fail


def main():
    ap = argparse.ArgumentParser(description="Batch submit Sfold structure prediction tasks")
    ap.add_argument("--type", choices=["WT", "MT", "BOTH"], default="BOTH",
                    help="Submit only WT / MT, or both (default)")
    ap.add_argument("--dry-run", action="store_true", help="Only list tasks to be submitted, do not open browser")
    args = ap.parse_args()

    types = ["WT", "MT"] if args.type == "BOTH" else [args.type]
    total_ok = total_fail = 0
    for stype in types:
        ok, fail = run_type(stype, dry_run=args.dry_run)
        total_ok += ok
        total_fail += fail
    if not args.dry_run:
        print("\n===== Summary =====  Success %d  Failed %d" % (total_ok, total_fail))


if __name__ == "__main__":
    main()
