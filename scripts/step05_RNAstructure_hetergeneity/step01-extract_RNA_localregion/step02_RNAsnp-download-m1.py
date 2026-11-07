# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
"""Batch-download RNAsnp Mode 1 result files from recorded URLs."""

import os
import traceback

import requests
from selenium import webdriver
from selenium.webdriver.common.by import By

element = "Intron"   # Element name
seqlen = "low"

alldir = "<LOCAL_DATA_ROOT>/negative-RNAsnp-rerun/Supple-RNAsnp-1/"
result_file = alldir + seqlen + "1wan_" + element + "_RNAsnp_result.txt"
out_base = "<LOCAL_DATA_ROOT>/negative-RNAsnp-rerun/"+element + "-" + seqlen + "-RNAsnp/"

# XPath, subdirectory, and extension for the 7 download targets on the result page
DOWNLOAD_TARGETS = [
    ("interval-p", ".txt",
     '//*[@id="description"]/table[2]/tbody/tr/td/a'),
    ("wt-mt-prob-pdf", ".pdf",
     '//*[@id="details"]/form/fieldset/table[2]/tbody/tr[3]/td/a[2]'),
    ("wt-mt-prob-value", ".txt",
     '//*[@id="details"]/form/fieldset/table[2]/tbody/tr[3]/td/a[3]'),
    ("wt-structure-pdf", ".pdf",
     '//*[@id="details"]/form/fieldset/table[4]/tbody/tr[2]/td/a[3]'),
    ("wt-structure-dp", ".txt",
     '//*[@id="details"]/form/fieldset/table[4]/tbody/tr[2]/td/a[4]'),
    ("mt-structure-pdf", ".pdf",
     '//*[@id="details"]/form/fieldset/table[4]/tbody/tr[5]/td/a[3]'),
    ("mt-structure-dp", ".txt",
     '//*[@id="details"]/form/fieldset/table[4]/tbody/tr[5]/td/a[4]'),
]


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


def download_all(driver, session, sample, url):
    """Visit one result page and download all 7 file types. Return whether all succeeded."""
    driver.get(url)

    for sub, ext, xpath in DOWNLOAD_TARGETS:
        out_dir = out_base + sub + "/"
        os.makedirs(out_dir, exist_ok=True)
        out_path = out_dir + sample + ext

        if os.path.isfile(out_path):
            print(f"  [Skip] {out_path} already exists")
            continue

        link = driver.find_element(By.XPATH, xpath).get_attribute("href")
        resp = session.get(link)
        with open(out_path, "wb") as f:
            f.write(resp.content)
        print(f"  [Downloaded] {out_path}")


def main():
    if not os.path.isfile(result_file):
        raise SystemExit(f"Result file does not exist: {result_file}")

    records = load_result_lines()
    print(f"{len(records)} records to process")

    # Reuse a single browser instance throughout
    driver = webdriver.Chrome()
    session = requests.Session()

    try:
        for sample, url in records:
            print(f"{sample}\t{url}")
            try:
                download_all(driver, session, sample, url)
            except Exception:
                print(f"[Failed] {sample} download exception, see error below")
                traceback.print_exc()
    finally:
        driver.quit()


if __name__ == "__main__":
    main()
