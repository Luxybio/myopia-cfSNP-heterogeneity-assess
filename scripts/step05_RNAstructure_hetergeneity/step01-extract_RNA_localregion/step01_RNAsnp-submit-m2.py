# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
"""Batch-submit local sequences to the RNAsnp web server (Mode 2)."""

import os
import time
import traceback

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

seqlen = "over"
element = "Intron"

BASE = "<LOCAL_DATA_ROOT>/negative-RNAsnp-rerun/Supple-RNAsnp-1/"
seq_dir = f"{BASE}/{element}-seq/{seqlen}1wan/"  # {m}/
snp_dir = f"{BASE}/{element}-snp/{seqlen}1wan/"
out_dir = f"{BASE}/{seqlen}1wan_{element}_RNAsnp_result.txt"  # {j}_

RNAsnp_URL = "https://rth.dk/resources/rnasnp/"
# Same as m1: manually specify local chromedriver (version matches Chrome), skip Selenium Manager
CHROMEDRIVER_PATH = (
    "<LOCAL_DATA_ROOT>/20261023-structure-paper-fighting/"
    "myopia_SNP_pyprogram/chromedriver134/chromedriver-win64/chromedriver.exe"
)
service = Service(executable_path=CHROMEDRIVER_PATH)


def load_done_records():
    """Load already submitted records for resume support. Returns a set of processed sample names."""
    done = set()
    if os.path.exists(out_dir):
        with open(out_dir, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    done.add(parts[1])
    return done


def append_result(idx, sample, url):
    """Append one result line. Retry on write failure; fall back to a backup file; finally print to avoid losing URLs."""
    line = f"{idx}\t{sample}\t{url}\n"
    for target in (out_dir, out_dir + ".fallback"):
        for attempt in range(3):
            try:
                with open(target, "a", encoding="utf-8") as f:
                    f.write(line)
                if target != out_dir:
                    print(f"[Info] Main result file write failed, written to fallback: {target}")
                return
            except OSError as e:
                print(f"[Warning] Failed to write {target} ({e}), attempt {attempt + 1}/3")
                time.sleep(1)
    print(f"[Result not persisted] {line.strip()}")


def submit_one(driver, rna_file, snp_file, jobname):
    """Fill and submit one task on the already opened RNAsnp page; return the result page URL."""
    driver.get(RNAsnp_URL)

    # Upload RNA sequence file and SNP file (verified form field names from the original script)
    driver.find_element(By.NAME, "RNASNP_a").send_keys(rna_file)
    driver.find_element(By.NAME, "RNASNP_t").send_keys(snp_file)

    # Select Mode 2 (local folding, RNAplfold): locate by name+value to avoid fragile long XPath
    driver.find_element(
        By.CSS_SELECTOR, "input[name='RNASNP_mode'][value='2']"
    ).click()

    # Job name
    driver.find_element(By.NAME, "job_name").send_keys(jobname)

    # Record URL before submission, then wait for the redirect to the result page
    before_url = driver.current_url
    driver.find_element(By.XPATH, "//button[normalize-space()='Submit']").click()

    # Explicit wait: up to 60 seconds for the URL to change (enter result page)
    WebDriverWait(driver, 60).until(EC.url_changes(before_url))
    time.sleep(2)  # Wait a bit more to ensure the result page URL is finalized
    return driver.current_url


def main():
    if not os.path.isdir(seq_dir):
        raise SystemExit(f"Sequence directory does not exist: {seq_dir}")

    done = load_done_records()
    files = sorted(f for f in os.listdir(seq_dir) if f.endswith(".txt"))
    print(f"Found {len(files)} sequence files, {len(done)} completed, "
          f"{len([f for f in files if os.path.splitext(f)[0] not in done])} remaining")

    # Reuse a single browser instance; start/close only once for the whole run
    # Same call as m1: webdriver.Chrome(service=service), confirmed to open pages correctly
    driver = webdriver.Chrome(service=service)
    driver.set_window_size(1200, 900)

    try:
        for idx, file in enumerate(files, start=1):
            sample = os.path.splitext(file)[0]
            if sample in done:
                print(f"[Skip] {sample} already in result file")
                continue

            rna_file = seq_dir + file
            snp_file = snp_dir + sample + ".txt"
            if not os.path.isfile(snp_file):
                print(f"[Warning] Missing corresponding SNP file, skip: {snp_file}")
                append_result("NA", sample, "MISSING_SNP_FILE")
                continue

            jobname = f"{idx}_{element}_{sample}"
            try:
                url = submit_one(driver, rna_file, snp_file, jobname)
                # Print first then persist: even if file write fails, the URL is not lost
                print(f"{idx}\t{sample}\t{url}", flush=True)
                append_result(idx, sample, url)
            except Exception:
                print(f"[Failed] {sample} submission exception, see error below")
                traceback.print_exc()
                append_result(idx, sample, "ERROR")
            finally:
                time.sleep(3)  # Submission interval to avoid server rate limits
    finally:
        driver.quit()


if __name__ == "__main__":
    main()
