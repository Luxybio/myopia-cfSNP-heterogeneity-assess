# myopia-cfSNP-heterogeneity-assess

Assess the heterogeneous effects of myopia GWAS susceptibility SNPs on genomic regulatory elements (GREs) and their RNA secondary structures.

## Overview

This repository contains a five-step analysis pipeline:

1. **SNP QC**
2. **SNP and GRE sequence extraction**
3. **Map SNPs onto GREs**
4. **TF heterogeneity**
5. **RNA structure heterogeneity**

## Repository Structure

```
.
├── data/              # example data and processed results
├── scripts/           # analysis scripts organized by step
├── requirements.txt   # Python dependencies
└── README.md
```

## Scripts

Scripts are organized by step under `scripts/`:

- `step01_SNPQC/` — SNP quality control
- `step02_SNPandGREsequence/` — sequence extraction
- `step03_SNPmappedonGRE/` — SNP mapping onto GREs
- `step04_TF_heterogeneity/` — TF binding analysis
- `step05_RNAstructure_hetergeneity/` — RNA structure analysis

## Setup

1. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. Replace `<LOCAL_DATA_ROOT>` and `<SERVER_DATA_ROOT>` in the scripts with your own paths.
3. Run the scripts in order from `step01` to `step05`.

## Environment

- **Python 3.x**: pandas, selenium, requests
- **R**: genetics, ggpubr, ggbreak, ggrepel, RNAsmc, RRNA
- **Tools**: bowtie2, Perl
- **Services**: RNAsnp, Sfold
