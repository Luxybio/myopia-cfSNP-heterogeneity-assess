# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
"""Merge RNAsnp interval p-value tables into the master element table."""

import csv
import re
from pathlib import Path

BASE = Path(r"<LOCAL_DATA_ROOT>")
NEG = BASE / "negative-RNAsnp-rerun"
element = "Exon"

ELEMENTS = {
    element: {
        "items": BASE / f"all-{element.lower()}-items.txt",
        "low_dir": NEG / f"{element}-low-RNAsnp" / "interval-p",
        "over_dir": NEG / f"{element}-over-RNAsnp" / "interval-p",
        "out": NEG / f"all{element.lower()}_p_interval.txt",
    }
}
    # Extensible (example):
    # "Exon": {
    #     "items": BASE / "all-exon-items.txt",
    #     "low_dir": NEG / "Exon-low-RNAsnp" / "interval-p",
    #     "over_dir": NEG / "Exon-over-RNAsnp" / "interval-p",
    #     "out": NEG / "allexon_p_interval.txt",
    # },


# 5 columns of the interval-p table (appended columns have the same names)
FIELDS = ["SNP", "FoldingWindow", "LocalRegion", "distance", "p-value"]


def read_items(path):
    """Read the header-containing master table (tab-delimited) and return (fieldnames, rows)."""
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        lines = [row for row in reader if any(c.strip() for c in row)]
    if not lines:
        raise SystemExit(f"Master table is empty: {path}")
    fieldnames = [c.strip() for c in lines[0]]
    # If duplicate column names appear in the header (e.g., two SNPonGRE columns),
    # add numeric suffixes to avoid dictionary key overwrites and lost values.
    seen = {}
    for i, c in enumerate(fieldnames):
        if c in seen:
            seen[c] += 1
            fieldnames[i] = f"{c}_{seen[c]}"
        else:
            seen[c] = 1
    if "SNPs" not in fieldnames or "SNPs_on_Genomic_regulatory_elements_fine" not in fieldnames:
        raise SystemExit(f"Master table missing SNPs / SNPs_on_Genomic_regulatory_elements_fine columns: {path}")
    rows = []
    for vals in lines[1:]:
        vals += [""] * (len(fieldnames) - len(vals))
        rows.append({k: (v or "").strip() for k, v in zip(fieldnames, vals)})
    return fieldnames, rows

# Filename: <Chr>_<rsID>_<allele>_<fine>.txt
FNAME_RE = re.compile(r"^(?P<chr>[^_]+)_(?P<snp>rs\d+)_(?P<allele>[AUCGTaucgt]+)_(?P<fine>.+)\.txt$")


def parse_filename(name):
    """Parse rsID and SNPs_on_Genomic_regulatory_elements_fine from the filename."""
    m = FNAME_RE.match(name)
    if not m:
        return None
    return m.group("snp"), m.group("fine")


def is_zip(path):
    """Check whether the file is a ZIP bundle (RNAsnp graphics bundles start with the PK header)."""
    try:
        with open(path, "rb") as fh:
            return fh.read(4) == b"PK\x03\x04"
    except OSError:
        return False


def extract_rows(path):
    """Read an interval-p file and return [(SNP, FoldingWindow, LocalRegion, distance, p_value), ...]."""
    records = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line[0] not in "AUCGaucg":   # Keep only data rows starting with A/U/C/G
                continue
            # A single line may contain multiple records (some result tables squeeze multiple rows into one line)
            for chunk in line.split(";"):
                chunk = chunk.strip()
                if not chunk:
                    continue
                cols = [c.strip() for c in chunk.split(",")]
                if len(cols) < 5:
                    print(f"  [Warning] Too few columns, skipping record: {chunk}")
                    continue
                records.append(tuple(cols[:5]))  # SNP, FoldingWindow, LocalRegion, distance, p-value
    return records


def process(element, cfg):
    items_path, out_path = cfg["items"], cfg["out"]
    if not items_path.is_file():
        raise SystemExit(f"[{element}] Cannot find master table: {items_path}")

    fieldnames, rows = read_items(items_path)
    if not fieldnames:
        raise SystemExit(f"[{element}] Master table is empty or missing header")

    # Build (SNPs, SNPs_on_Genomic_regulatory_elements_fine) -> list of row indices
    index = {}
    for i, row in enumerate(rows):
        k = (row.get("SNPs", "").strip(),
             row.get("SNPs_on_Genomic_regulatory_elements_fine", "").strip())
        index.setdefault(k, []).append(i)

    new_cols = FIELDS[:]   # low/over not distinguished; column names match interval-p table
    hits = {}              # row index -> record list (low/over merged)
    unmatched_files, zip_files = [], []

    for tag, d in (("low", cfg["low_dir"]), ("over", cfg["over_dir"])):
        if not d.is_dir():
            print(f"[{element}] [Warning] Directory not found: {d}, {tag} results will be left blank")
            continue
        files = sorted(p for p in d.iterdir() if p.suffix.lower() == ".txt")
        print(f"[{element}] {tag}: found {len(files)} interval-p files")

        for path in files:
            if is_zip(path):
                print(f"[Skip] ZIP graphics bundle (no p-value table): {path.name}")
                zip_files.append(path.name)
                continue

            parsed = parse_filename(path.name)
            if parsed is None:
                print(f"[Warning] Filename does not match expected pattern, skipped: {path.name}")
                unmatched_files.append(path.name)
                continue
            snp_id, fine = parsed
            key = (snp_id, fine)
            if key not in index:
                print(f"[Warning] No matching row in master table: {path.name}  ->  {key}")
                unmatched_files.append(path.name)
                continue

            records = extract_rows(path)
            if not records:
                print(f"[Warning] No A/U/C/G-starting data rows extracted: {path.name}")
                unmatched_files.append(path.name)
                continue

            for i in index[key]:
                hits.setdefault(i, []).extend(records)

    # Write output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames + new_cols, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for i, row in enumerate(rows):
            out = dict(row)
            recs = hits.get(i)
            for c in FIELDS:
                out[c] = ";".join(r[FIELDS.index(c)] for r in recs) if recs else ""
            writer.writerow(out)

    filled = len(hits)
    print(f"\n[{element}] Output file: {out_path}")
    print(f"[{element}] Master table rows: {len(rows)} | Rows filled: {filled} "
          f"| Unmatched files: {len(unmatched_files)} | ZIP graphics bundles (no p table): {len(zip_files)}")
    if unmatched_files:
        print("Unmatched file list:")
        for n in unmatched_files:
            print("  -", n)
    if zip_files:
        print(f"Tip: {len(zip_files)} files are ZIP graphics bundles (RNAsnp only saved graphic output); "
              f"resubmit those jobs to obtain p-value tables and rerun this script.")


def main():
    for element, cfg in ELEMENTS.items():
        print(f"\n========== Processing element: {element} ==========")
        process(element, cfg)


if __name__ == "__main__":
    main()
