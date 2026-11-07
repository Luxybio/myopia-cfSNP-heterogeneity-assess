# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
# RNAsnp web server batch submission script (Selenium 4.x compatible)
# Fixes:
#   1) webdriver.Chrome no longer accepts a chromedriver path argument; Selenium Manager auto-matches the driver
#   2) find_element_by_name -> find_element(By.NAME, ...)
#   3) Initialize i
import os
import time

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By

# Manually specify chromedriver (downloaded into the project directory, version 134.0.6998.165 matches Chrome)
CHROMEDRIVER_PATH = (
    "<LOCAL_DATA_ROOT>/20261023-structure-paper-fighting/"
    "myopia_SNP_pyprogram/chromedriver134/chromedriver-win64/chromedriver.exe"
)
service = Service(executable_path=CHROMEDRIVER_PATH)

seqlen = "low"
element = "Intron"
dir = "<LOCAL_DATA_ROOT>/negative-RNAsnp-rerun/Supple-RNAsnp-1/" + element + "-seq/" + seqlen + "1wan/"
# SNP directory structure mirrors sequence directory, with low1wan/over1wan subdirectories below
snpdir = "<LOCAL_DATA_ROOT>/negative-RNAsnp-rerun/Supple-RNAsnp-1/" + element + "-snp/" + seqlen + "1wan/"
outdir = "<LOCAL_DATA_ROOT>/negative-RNAsnp-rerun/Supple-RNAsnp-1/" + seqlen + "1wan_" + element + "_RNAsnp_result.txt"

filedir = os.listdir(dir)

i = 0  # Initialize counter, otherwise i=i+1 raises NameError
for file in filedir:
    if not file.endswith(".txt"):
        continue

    curfile = file.split(".txt")[0]
    rnafile = dir + curfile + ".txt"
    snpfile = snpdir + curfile + ".txt"

    # Manually specify chromedriver path; do not rely on automatic download
    driver = webdriver.Chrome(service=service)
    try:
        driver.get("https://rth.dk/resources/rnasnp/")

        # Upload sequence file (input type=file; send_keys the absolute path)
        element_rna = driver.find_element(By.NAME, "RNASNP_a")
        element_rna.send_keys(rnafile)

        # Upload SNP file
        element_snp = driver.find_element(By.NAME, "RNASNP_t")
        element_snp.send_keys(snpfile)

        # Job name
        i = i + 1
        jobname = str(i) + "_" + element + "_" + curfile
        element_jobname = driver.find_element(By.NAME, "job_name")
        element_jobname.send_keys(jobname)

        # Submit button (page shows <button onclick="submit1()">Submit</button>)
        element_search_button = driver.find_element(By.CSS_SELECTOR, "button[onclick]")
        element_search_button.click()

        time.sleep(20)  # Wait for redirect to result page
        url = driver.current_url
    finally:
        driver.quit()  # Ensure browser closes even on error to avoid leftover processes

    output = open(outdir, "a")
    output.write(str(i) + "\t" + curfile + "\t" + url + "\n")
    output.close()
    print(str(i) + "\t" + curfile + "\t" + url + "\n")
