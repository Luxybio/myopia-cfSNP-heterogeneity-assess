#!/usr/bin/env python3
# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
"""Extract WT/MT local sequences for Sfold structure prediction."""

import argparse
import os
import re
import sys
from collections import defaultdict

# Default paths (modify only here; no need to touch logic below)
DEFAULT_ITEMS     = r"<LOCAL_DATA_ROOT>\Sfold-seq-build\Sfold-rerun-items.txt"
# Multiple sequence library root directories are supported (searched in order, first match wins)
DEFAULT_SEQ_ROOTS = [
    r"<LOCAL_DATA_ROOT>\negative-RNAsnp-rerun",
    r"<LOCAL_DATA_ROOT>\negtive-RNAsnp-16-rerun",   # Other half for Exon-seq etc.
]
DEFAULT_WT_DIR    = r"<LOCAL_DATA_ROOT>\Sfold-seq-build\WT"
DEFAULT_MT_DIR    = r"<LOCAL_DATA_ROOT>\Sfold-seq-build\MT"

SEQ_DIR_GLOB = "*-seq"          # Subdirectory wildcard participating in the search
SEQ_EXT      = ".txt"

# Column names (support smart distinction of duplicate column names)
COL_RSID   = "SNPs"
COL_GREID  = "SNPs_on_Genomic_regulatory_elements"
COL_RANGE  = "range"
COL_SNP    = "SNP"
COL_REF    = "Ref_allele"
COL_ALT    = "Alt_allele"
COL_REGELE = "Regulatory_elements"
COL_CHR    = "Chr"
COL_LEN    = "length"
COL_GENE   = "Gene_symbol"

COMPLEMENT = str.maketrans(
    "ACGTUNRYKMacgtunrykm",
    "TGCAANYRMKtgcaanyrmk",
)

LOG = []          # Per-entry results, written to log at the end


# Basic utilities
def norm_allele(s):
    """Normalize allele: strip whitespace, uppercase, U -> T (original sequence files use T)."""
    return (s or "").strip().upper().replace("U", "T")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    """Read single/multiple FASTA entries, returning [(header, seq), ...]; also tolerate bare sequence files without '>'."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        raw = fh.read()
    records, header, buf = [], None, []
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(buf)))
            header, buf = line[1:].strip(), []
        else:
            if header is None:
                header = os.path.basename(path)
            buf.append(line)
    if header is not None:
        records.append((header, "".join(buf)))
    # Keep only segments containing ACGTUN, removing possible spaces/digits
    clean = []
    for h, s in records:
        s = re.sub(r"[^ACGTUNacgtun]", "", s).upper()
        if s:
            clean.append((h, s))
    return clean


# Items table parsing
def split_header(header):
    """Split header into [(colname, colindex), ...], keeping all positions for duplicate names."""
    return [(name.strip(), i) for i, name in enumerate(header)]


def pick_columns(cols):
    """
    Resolve the column indices needed.
    SNPs_on_Genomic_regulatory_elements has two columns:
      - pure digits = SNP position within this GRE sequence;
      - contains ENST/ENSE = GRE ID used as the search keyword.
    """
    idx = {}
    idx["_gre_candidates"] = [i for name, i in cols if name == COL_GREID]

    for key, name in (("rsid", COL_RSID), ("range", COL_RANGE), ("snp", COL_SNP),
                      ("ref", COL_REF), ("alt", COL_ALT), ("elem", COL_REGELE),
                      ("chr", COL_CHR), ("len", COL_LEN), ("gene", COL_GENE)):
        hit = [i for n, i in cols if n == name]
        idx[key] = hit[0] if hit else None
    return idx


def classify_gre_columns(rows, candidates):
    """Distinguish GRE ID columns from numeric position columns by content."""
    gre_id_idx, gre_pos_idx = None, None
    if not candidates:
        return None, None
    if len(candidates) == 1:
        sample = rows[0][candidates[0]] if rows else ""
        if re.fullmatch(r"\d+", sample.strip()):
            return None, candidates[0]
        return candidates[0], None
    # Multiple columns: the one containing letters is the ID column
    for i in candidates:
        vals = [r[i].strip() for r in rows if i < len(r) and r[i].strip()]
        if vals and sum(1 for v in vals if re.search(r"[A-Za-z]", v)) > len(vals) / 2:
            gre_id_idx = i
        else:
            gre_pos_idx = i
    return gre_id_idx, gre_pos_idx


def parse_range(txt):
    """'29484-29552' -> (29484, 29552), 1-based inclusive."""
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", txt or "")
    if not m:
        return None, None
    a, b = int(m.group(1)), int(m.group(2))
    return (a, b) if a <= b else (b, a)


def parse_snp(txt):
    """
    'G29504C' -> ('G', 29504, 'C')    wild-type allele / 1-based position / mutant allele
    Also supports 'c.29504G>C', '29504 G C', etc.
    """
    if not txt:
        return None
    t = txt.strip().replace(" ", "")
    m = re.fullmatch(r"([A-Za-z]+)(\d+)([A-Za-z]+)", t)
    if m:
        return norm_allele(m.group(1)), int(m.group(2)), norm_allele(m.group(3))
    m = re.search(r"(\d+)\s*([A-Za-z])>([A-Za-z])", t)
    if m:
        return norm_allele(m.group(3)), int(m.group(1)), norm_allele(m.group(2))
    m = re.search(r"(\d+)\s*([ACGTUacgtu])\s*([ACGTUacgtu])", t)
    if m:
        return norm_allele(m.group(2)), int(m.group(1)), norm_allele(m.group(3))
    return None


# Sequence file indexing and search
def build_seq_index(seq_roots):
    """
    Recursively find all directories matching *-seq under the given seq_roots (any depth,
    e.g., xxx/Supple-RNAsnp/Intron-seq/low1wan) and collect all *.txt sequence files.
    Returns {basename_without_ext: [fullpath, ...]}
    """
    index = defaultdict(list)
    seq_re = re.compile(SEQ_DIR_GLOB.replace("*", ".*"))
    for seq_root in seq_roots:
        if not os.path.isdir(seq_root):
            sys.exit("[FATAL] Sequence root directory does not exist: %s" % seq_root)
        for dirpath, dirnames, filenames in os.walk(seq_root):
            if not seq_re.fullmatch(os.path.basename(dirpath)):
                continue
            # Hit *-seq directory: collect files from it and all its descendants (low1wan/over1wan etc.)
            for sub, _sd, fs in os.walk(dirpath):
                for fn in fs:
                    if not fn.lower().endswith(SEQ_EXT):
                        continue
                    stem = fn[: -len(SEQ_EXT)]
                    index[stem].append(os.path.join(sub, fn))
            dirnames[:] = []   # Already traversed deeply above; do not descend again in outer loop
    return index


def parse_seq_filename_basename(stem):
    """Filename {Chr}_{rsID}_{RefAllele}_{GREid}.txt -> (chr, rsid, ref, greid)"""
    parts = stem.split("_", 3)
    if len(parts) < 4:
        return None
    return parts[0], parts[1], parts[2], parts[3]


def candidate_score(stem, path, rsid, ref, gre_id, decl_len, snp_pos, wt):
    """
    Score a candidate file:
      +6  GRE ID is an exact match (e.g., gid has no segment suffix and file also has none)
      +3  Sequence length matches the declared length in items
      +5  SNP base equals wild-type allele
      +2  SNP base equals the complement of wild-type allele (sequence may be reverse complemented)
      +1  GRE gene ID matches (only segment suffix differs)
      -90 SNP base is neither wild-type nor its complement
    """
    info = parse_seq_filename_basename(stem)
    if info is None:
        return None, None
    f_chr, f_rsid, f_ref, f_gid = info
    if f_rsid != rsid:
        return None, None
    if f_gid.split("-")[0] != gre_id.split("-")[0]:
        return None, None

    score = 0
    if f_gid == gre_id:
        score += 6
    recs = read_fasta(path)
    if not recs:
        return None, None
    seq = recs[0][1]
    if decl_len and len(seq) == decl_len:
        score += 3
    base = seq[snp_pos - 1] if 0 < snp_pos <= len(seq) else ""
    if base and base == wt:
        score += 5
    elif base and base == wt.translate(COMPLEMENT):
        score += 2
    elif base:
        score -= 90
    return score, (stem, path, seq, base, recs[0][0])


def locate_seq(index, rsid, gre_id, ref, decl_len, snp_pos, wt):
    """
    Keyword search: rsID + GRE ID (ENST/ENSE).
    Prefer exact-name files, otherwise choose by score.
    """
    # 1) Exact GRE ID match (e.g., ENST00000640003-22765277-22822017, ENST00000486151)
    exact = []
    for stem, paths in index.items():
        info = parse_seq_filename_basename(stem)
        if not info:
            continue
        _c, f_rsid, _r, f_gid = info
        if f_rsid == rsid and f_gid == gre_id:
            for p in paths:
                exact.append((stem, p))
    pool = exact
    if not pool:  # 2) Fall back to same gene ID prefix
        for stem, paths in index.items():
            info = parse_seq_filename_basename(stem)
            if not info:
                continue
            _c, f_rsid, _r, f_gid = info
            if f_rsid == rsid and f_gid.split("-")[0] == gre_id.split("-")[0]:
                for p in paths:
                    pool.append((stem, p))

    scored = []
    for stem, p in pool:
        s, payload = candidate_score(stem, p, rsid, ref, gre_id, decl_len, snp_pos, wt)
        if payload is not None:
            scored.append((s, payload))
    if not scored:
        return None
    scored.sort(key=lambda x: (-x[0], len(x[1][2])))
    return scored[0][1]


# Main workflow
def load_items(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = [l.rstrip("\n").rstrip("\r") for l in fh if l.strip("\r\n")]
    if not lines:
        sys.exit("[FATAL] items file is empty: %s" % path)
    header = lines[0].split("\t")
    cols = split_header(header)
    idx = pick_columns(cols)
    rows = [l.split("\t") for l in lines[1:]]
    gre_id_idx, gre_pos_idx = classify_gre_columns(rows, idx["_gre_candidates"])
    idx["greid"], idx["grepos"] = gre_id_idx, gre_pos_idx
    return idx, rows


def write_fasta(path, header, seq):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(">%s\n" % header)
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")


def main():
    ap = argparse.ArgumentParser(description="Extract Sfold local WT/MT sequences")
    ap.add_argument("--items", default=DEFAULT_ITEMS)
    ap.add_argument("--seq-root", nargs="+", default=DEFAULT_SEQ_ROOTS,
                    help="Sequence library root directories; multiple roots can be provided and will be searched in order")
    ap.add_argument("--wt", default=DEFAULT_WT_DIR)
    ap.add_argument("--mt", default=DEFAULT_MT_DIR)
    ap.add_argument("--dry-run", action="store_true", help="Check only, do not write files")
    ap.add_argument("--on-outside", choices=["keep", "shift", "skip"], default="keep",
                    help="How to handle SNP positions falling outside the range: "
                         "keep=keep original range (MT same as WT, default) / "
                         "shift=keep window length and center window on the SNP / "
                         "skip=skip this entry")
    args = ap.parse_args()
    opts = args

    idx, rows = load_items(args.items)
    if idx["greid"] is None:
        sys.exit("[FATAL] Could not locate GRE ID column (%s) in header" % COL_GREID)
    print("[INFO] items rows: %d" % len(rows))
    print("[INFO] GRE ID column = column %d, GRE position column = %s"
          % (idx["greid"], idx["grepos"]))

    index = build_seq_index(args.seq_root)
    print("[INFO] Sequence index built: %d unique file names" % len(index))

    if not args.dry_run:
        os.makedirs(args.wt, exist_ok=True)
        os.makedirs(args.mt, exist_ok=True)

    ok, skip, warn = 0, 0, 0
    used_names = defaultdict(int)

    for n, r in enumerate(rows, start=1):
        def col(key):
            i = idx.get(key)
            return r[i].strip() if (i is not None and i < len(r)) else ""

        rsid = col("rsid")
        gre_id = r[idx["greid"]].strip() if idx["greid"] < len(r) else ""
        rng_txt = col("range")
        snp_txt = col("snp")
        ref = col("ref")
        alt = col("alt")
        elem = col("elem")
        chr_ = col("chr")
        gene = col("gene")
        decl_len_txt = col("len")
        decl_len = int(decl_len_txt) if decl_len_txt.isdigit() else 0

        tag = "[%03d] Chr%s %s %s %s" % (n, chr_, rsid, elem or "-", gre_id)

        # ---- Basic validation ----
        if not rsid or not gre_id:
            LOG.append("%s\tSKIP\tmissing rsID or GRE ID" % tag)
            skip += 1
            continue

        start, end = parse_range(rng_txt)
        if start is None:
            LOG.append("%s\tSKIP\trange column cannot be parsed: %r" % (tag, rng_txt))
            skip += 1
            continue

        snp = parse_snp(snp_txt)
        if snp is None:
            LOG.append("%s\tSKIP\tSNP column cannot be parsed: %r" % (tag, snp_txt))
            skip += 1
            continue
        wt, snp_pos, mt = snp          # wt/mt already U->T; snp_pos is 1-based within the GRE sequence
        # The SNP column gives the allele on the transcript strand; Ref_allele is on the genomic plus strand.
        # For negative-strand genes the two are complements, which is normal;
        # only flag when they are neither identical nor complementary.
        ref_n = norm_allele(ref)
        if ref_n and ref_n != wt and ref_n != wt.translate(COMPLEMENT):
            LOG.append("%s\tWARN\tSNP wild-type allele(%s) and Ref_allele(%s) are neither identical nor complementary; please review"
                       % (tag, wt, ref_n))
            warn += 1

        # ---- Find original sequence file ----
        hit = locate_seq(index, rsid, gre_id, wt, decl_len, snp_pos, wt)
        if hit is None:
            LOG.append("%s\tSKIP\tNo matching -seq sequence file found (keywords: %s + %s)"
                       % (tag, rsid, gre_id))
            skip += 1
            continue
        stem, src_path, full_seq, base, src_header = hit

        # ---- Out-of-bounds check ----
        L = len(full_seq)
        if end > L:
            LOG.append("%s\tWARN\trange end %d exceeds sequence length %d, truncated" % (tag, end, L))
            warn += 1
            end = L
        if start > L:
            LOG.append("%s\tSKIP\trange start %d exceeds sequence length %d" % (tag, start, L))
            skip += 1
            continue

        # ---- Determine sequence orientation (transcript strand or reverse complement) ----
        strand_flip = False
        if base == wt:
            pass
        elif base and base == wt.translate(COMPLEMENT):
            strand_flip = True
            LOG.append("%s\tWARN\tSite base is the complement of wild-type(%s/%s); treat as reverse complement" % (tag, base, wt))
            warn += 1
        else:
            LOG.append("%s\tWARN\tSite base(%s) does not match wild-type(%s); still replacing by wild-type, please review manually"
                       % (tag, base or "?", wt))
            warn += 1

        # ---- Whether SNP falls inside the window; apply --on-outside if not ----
        if not (start <= snp_pos <= end):
            win_len = end - start + 1
            msg = "%s\tWARN\tSNP position %d is not within range %d-%d" % (tag, snp_pos, start, end)
            if opts.on_outside == "skip":
                LOG.append(msg + ", skipped (--on-outside skip)")
                skip += 1
                continue
            if opts.on_outside == "shift":
                new_start = max(1, snp_pos - win_len // 2)
                new_end = min(L, new_start + win_len - 1)
                new_start = max(1, new_end - win_len + 1)
                LOG.append(msg + ", window shifted to %d-%d (--on-outside shift)"
                           % (new_start, new_end))
                start, end = new_start, new_end
            else:
                LOG.append(msg + ", MT same as WT (--on-outside keep)")
            warn += 1

        # ---- Extract local sequence ----
        wt_local = full_seq[start - 1:end]
        mt_base = mt.translate(COMPLEMENT) if strand_flip else mt
        off = snp_pos - start          # 0-based offset
        if 0 <= off < len(wt_local):
            mt_local = wt_local[:off] + mt_base + wt_local[off + 1:]
            if wt_local[off] != (wt.translate(COMPLEMENT) if strand_flip else wt):
                LOG.append("%s\tWARN\tLocal sequence base at position %d (%s) does not match expected wild-type"
                           % (tag, off + 1, wt_local[off]))
                warn += 1
        else:
            mt_local = wt_local
            if opts.on_outside == "shift":
                LOG.append("%s\tWARN\tAfter shifting, SNP position %d still not within range %d-%d"
                           % (tag, snp_pos, start, end))
            warn += 1

        # ---- Output (file names reuse original sequence file names) ----
        out_name = stem + SEQ_EXT
        if used_names[out_name]:
            new_name = "%s__dup%d%s" % (stem, used_names[out_name] + 1, SEQ_EXT)
            LOG.append("%s\tWARN\tDuplicate output file name, renamed to %s" % (tag, new_name))
            warn += 1
            out_name = new_name
        used_names[out_name] += 1

        header = src_header or "%s-%d-%d" % (gre_id, start, end)   # Reuse original file header
        wt_out = os.path.join(args.wt, out_name)
        mt_out = os.path.join(args.mt, out_name)
        if not args.dry_run:
            write_fasta(wt_out, header, wt_local)
            write_fasta(mt_out, header, mt_local)

        src_root = next((rt for rt in args.seq_root if src_path.startswith(os.path.abspath(rt))),
                        args.seq_root[0])
        LOG.append("%s\tOK\tsrc=%s\trange=%d-%d\tlen=%d\tSNP=%s%d%s(position %d in window)\tref/alt=%s/%s"
                   % (tag, os.path.relpath(src_path, src_root), start, end,
                      len(wt_local), wt, snp_pos, mt, off + 1 if 0 <= off < len(wt_local) else -1,
                      ref or "-", alt or "-"))
        ok += 1

    # ---- Write log ----
    log_path = os.path.join(os.path.dirname(os.path.abspath(args.items)),
                            "build-WT-MT-localseq.log.txt")
    if not args.dry_run:
        with open(log_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("build-WT-MT-localseq.py run log\n")
            fh.write("items: %s\nWT: %s\nMT: %s\n\n" % (args.items, args.wt, args.mt))
            fh.write("\n".join(LOG) + "\n")

    print("\n===== Summary =====")
    print("Success: %d    Skipped: %d    Warnings: %d    Total rows: %d" % (ok, skip, warn, len(rows)))
    if warn:
        print("\n--- Warning details ---")
        for l in LOG:
            if "\tWARN\t" in l or "\tNOTE\t" in l:
                print("  " + l)
    if skip:
        print("\n--- Skip details ---")
        for l in LOG:
            if "\tSKIP\t" in l:
                print("  " + l)
    if not args.dry_run:
        print("\n[INFO] Log: %s" % log_path)
        print("[INFO] WT  -> %s" % args.wt)
        print("[INFO] MT  -> %s" % args.mt)


if __name__ == "__main__":
    main()
