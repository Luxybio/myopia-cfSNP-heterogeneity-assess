#!/usr/bin/env python3
# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
"""Select the most likely Sfold centroid using RNAsnp pairing probabilities."""

import argparse
import glob
import os
import re
import shutil
import sys

# Default configuration
DEFAULT_DOWN_ROOT = r"<LOCAL_DATA_ROOT>\Sfold-seq-build\sfold-down"
DEFAULT_PROB_ROOTS = [
    r"<LOCAL_DATA_ROOT>\negative-RNAsnp-rerun",
    r"<LOCAL_DATA_ROOT>\negtive-RNAsnp-16-rerun",
]
DEFAULT_OUT_ROOT = r"<LOCAL_DATA_ROOT>\20261023-structure-paper-fighting\myopia_SNP_pyprogram\sfold-selectstructrue-4"

PROB_SUBDIR = "wt-mt-prob-value"      # Subdirectory containing RNAsnp probability matrices
SEQ_EXT = ".ct"
TYPES = ("WT", "MT")
PROB_TRIANGLE = {"WT": "upper", "MT": "lower"}   # See header comment

LOG = []


def log(msg):
    LOG.append(msg)
    print(msg, flush=True)


# Read Sfold .ct
def parse_ct(path):
    """
    Return dict:
        length / energy(mfe value) / init_energy / pairs(set of (i,j), i<j) / kind
    Sfold .ct columns: index base prev next pairing stack
    """
    length = energy = init_energy = None
    pairs = set()
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            if length is None and "Length" in line:
                m = re.search(r"Length:\s*(\d+)", line)
                if m:
                    length = int(m.group(1))
            if energy is None and "Energy" in line:
                m = re.search(r"Energy:\s*(-?\d+\.?\d*)", line)
                if m:
                    energy = float(m.group(1))
            m = re.search(r"initially\s*(-?\d+\.?\d*)", line)
            if init_energy is None and m:
                init_energy = float(m.group(1))
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                i = int(fields[0])
                j = int(fields[4])
            except ValueError:
                continue
            if j > 0 and j != i:
                pairs.add((min(i, j), max(i, j)))
    return {"length": length, "energy": energy, "init_energy": init_energy,
            "pairs": pairs, "name": os.path.basename(path), "path": path}


# Probability matrix index
def build_prob_index(prob_roots):
    """{jobkey: [matrix file paths, ...]}, keeping all duplicates (e.g., Supple copies)"""
    idx = {}
    for root in prob_roots:
        if not os.path.isdir(root):
            log("[WARN] Probability root directory does not exist: %s" % root)
            continue
        pat = os.path.join(root, "**", PROB_SUBDIR, "*.txt")
        for p in glob.glob(pat, recursive=True):
            key = os.path.basename(p)[: -len(".txt")]
            idx.setdefault(key, []).append(p)
    return idx


def load_matrix(path):
    """Read an N×N probability matrix (tab-delimited; trailing empty columns may exist)"""
    mat = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            row = []
            for x in line.split("\t"):
                x = x.strip()
                if x == "":
                    continue
                try:
                    row.append(float(x))
                except ValueError:
                    row.append(0.0)
            if row:
                mat.append(row)
    n = len(mat)
    for r in mat:                      # Align row lengths
        if len(r) < n:
            r.extend([0.0] * (n - len(r)))
    return [r[:n] for r in mat]


def pair_prob_lookup(mat, typ):
    """Return {(i,j): p} (i<j). WT uses the lower triangle, MT uses the upper triangle."""
    n = len(mat)
    lut = {}
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            try:
                v = mat[i - 1][j - 1] if PROB_TRIANGLE[typ] == "upper" else mat[j - 1][i - 1]
            except IndexError:
                v = 0.0
            if v > 0:
                lut[(i, j)] = v
    return lut


def pick_prob_file(paths, expect_len):
    """If multiple matrices files with the same name exist (original + Supple copies): prefer the one whose row count matches the Sfold window"""
    if not paths:
        return None
    for p in paths:
        try:
            n = sum(1 for _ in open(p, "r", encoding="utf-8", errors="replace"))
        except OSError:
            continue
        if n == expect_len:
            return p
    return paths[0]


# Single task processing
def list_centroids(struct_dir, include_ensemble=True):
    out = []
    for f in sorted(os.listdir(struct_dir)):
        low = f.lower()
        if not low.endswith(SEQ_EXT):
            continue
        if "ecentroid" in low:
            if include_ensemble:
                out.append(os.path.join(struct_dir, f))
            continue
        if "centroid" in low:          # c01.ccentroid.ct ...
            out.append(os.path.join(struct_dir, f))
    return out


def process_job(typ, jobname, job_dir, prob_paths, out_root,
                min_prob=0.0, include_ensemble=True, dry_run=False):
    """Return the result dict for this task, or None on failure."""
    jobkey = re.sub(r"^%s-" % typ, "", jobname)
    subdirs = [d for d in os.listdir(job_dir) if os.path.isdir(os.path.join(job_dir, d))]
    if not subdirs:
        log("[SKIP] %s %s: no <pid> subdirectory (download may be incomplete)" % (typ, jobname))
        return None
    pid = sorted(subdirs)[0]
    struct_dir = os.path.join(job_dir, pid, "structures")
    if not os.path.isdir(struct_dir):
        log("[SKIP] %s %s: no structures directory" % (typ, jobname))
        return None

    cts = list_centroids(struct_dir, include_ensemble=include_ensemble)
    if not cts:
        log("[SKIP] %s %s: no centroid structures under structures" % (typ, jobname))
        return None

    # Probability matrix
    lut = None
    prob_src = ""
    ct_len = parse_ct(cts[0])["length"]
    if prob_paths:
        cand = pick_prob_file(prob_paths, ct_len)
        prob_src = cand
        try:
            mat = load_matrix(cand)
            if ct_len and len(mat) != ct_len:
                log("[WARN] %s %s: probability matrix dimension %d ≠ Sfold window length %s, coordinates may be misaligned"
                    % (typ, jobname, len(mat), ct_len))
            lut = pair_prob_lookup(mat, typ)
        except Exception as exc:                                  # noqa: BLE001
            log("[WARN] %s %s: failed to read probability matrix (%s), falling back to lowest energy selection" % (typ, jobname, exc))
            lut = None
    if lut is None:
        log("[WARN] %s %s: no usable probability matrix, falling back to lowest ΔG selection" % (typ, jobname))

    rows = []
    for ct in cts:
        info = parse_ct(ct)
        pairs = info["pairs"]
        if lut is not None:
            hit = [lut[k] for k in pairs if k in lut and lut[k] > min_prob]
            n_hit, psum = len(hit), sum(hit)
        else:
            n_hit, psum = 0, 0.0
        rows.append({"ct": info["name"], "n_hit": n_hit, "psum": psum,
                     "npair": len(pairs), "length": info["length"], "pairs": pairs,
                     "energy": info["energy"], "init": info["init_energy"],
                     "path": ct})

    if lut is not None:
        rows.sort(key=lambda r: (-r["n_hit"], -r["psum"],
                                 r["energy"] if r["energy"] is not None else 0.0, r["ct"]))
    else:
        rows.sort(key=lambda r: (r["energy"] if r["energy"] is not None else 0.0, r["ct"]))

    best = rows[0]
    ct_stem = os.path.splitext(best["ct"])[0].replace(".", "-")
    out_name = "%s_%s_%s%s" % (typ, jobkey, ct_stem, SEQ_EXT)   # Naming follows old program
    out_dir = os.path.join(out_root, typ)
    out_path = os.path.join(out_dir, out_name)

    if not dry_run:
        os.makedirs(out_dir, exist_ok=True)
        shutil.copy(best["path"], out_path)
        # Candidate table per task (old selectstructure style)
        ss_dir = os.path.join(out_root, "selectstructure", typ)
        os.makedirs(ss_dir, exist_ok=True)
        ss_path = os.path.join(ss_dir, "%s_%s-selectstructure.txt" % (typ, jobkey))
        with open(ss_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("centroid\tmatched_pairs\tsum_prob\tpairs\tdG\tselected\n")
            for r in rows:
                fh.write("%s\t%d\t%.6f\t%d\t%s\t%s\n"
                         % (r["ct"], r["n_hit"], r["psum"], r["npair"],
                            "" if r["energy"] is None else "%.2f" % r["energy"],
                            "1" if r is best else "0"))

    log("[OK] %s %-60s pid=%s selected %s (hits %d/%d pairs, sum p=%.3f, ΔG=%s)"
        % (typ, jobname, pid, best["ct"], best["n_hit"], best["npair"], best["psum"],
           "NA" if best["energy"] is None else "%.2f" % best["energy"]))

    return {"type": typ, "job": jobname, "jobkey": jobkey, "pid": pid,
            "prob_src": prob_src, "rows": rows, "best": best,
            "out_name": out_name, "out_path": out_path}


# Main workflow
def run(down_root, prob_roots, out_root, types, min_prob, include_ensemble, dry_run):
    prob_idx = build_prob_index(prob_roots)
    log("[INFO] Probability matrix index: %d jobkeys (from %d root directories)"
        % (len(prob_idx), len(prob_roots)))

    results = []
    for typ in types:
        type_dir = os.path.join(down_root, typ)
        if not os.path.isdir(type_dir):
            log("[WARN] Directory does not exist: %s" % type_dir)
            continue
        jobs = sorted(d for d in os.listdir(type_dir)
                      if os.path.isdir(os.path.join(type_dir, d)))
        log("[INFO] %s: %d tasks" % (typ, len(jobs)))
        for n, jobname in enumerate(jobs, 1):
            jobkey = re.sub(r"^%s-" % typ, "", jobname)
            prob_paths = prob_idx.get(jobkey, [])
            if not prob_paths:
                log("[WARN] %s %s: no probability matrix with the same name (jobkey=%s)" % (typ, jobname, jobkey))
            print("  (%d/%d) %s" % (n, len(jobs), jobname), flush=True)
            res = process_job(typ, jobname, os.path.join(type_dir, jobname), prob_paths,
                              out_root, min_prob=min_prob,
                              include_ensemble=include_ensemble, dry_run=dry_run)
            if res:
                results.append(res)

    if dry_run:
        log("[DRY-RUN] No files written")
        return results

    os.makedirs(out_root, exist_ok=True)

    # 1) Aggregated long table of all centroids
    p1 = os.path.join(out_root, "centroid-prob-count.txt")
    with open(p1, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("type\tjobname\tpid\tcentroid\tmatched_pairs\tsum_prob\tpairs\tlength\tdG\tselected\n")
        for r in results:
            for row in r["rows"]:
                fh.write("%s\t%s\t%s\t%s\t%d\t%.6f\t%d\t%s\t%s\t%s\n"
                         % (r["type"], r["job"], r["pid"], row["ct"], row["n_hit"],
                            row["psum"], row["npair"], row["length"],
                            "" if row["energy"] is None else "%.2f" % row["energy"],
                            "1" if row is r["best"] else "0"))

    # 2) Selected structure energy table (equivalent to old program's centroid-energy.txt)
    p2 = os.path.join(out_root, "centroid-energy.txt")
    with open(p2, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("ct_file\tlength\tdG_kcal_per_mol\tMFE_of_input\tsource_centroid\tjobname\n")
        for r in results:
            b = r["best"]
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\n"
                     % (r["out_name"], b["length"],
                        "" if b["energy"] is None else "%.2f" % b["energy"],
                        "" if b["init"] is None else "%.2f" % b["init"],
                        b["ct"], r["job"]))

    # 3) WT / MT pairing comparison
    p3 = os.path.join(out_root, "WT-MT-centroid-compare.txt")
    index = {(r["type"], r["jobkey"]): r for r in results}
    with open(p3, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("jobkey\tWT_centroid\tWT_matched\tWT_dG\tMT_centroid\tMT_matched\tMT_dG\t"
                 "dDG_MT_minus_WT\tpair_changed\tpair_common\n")
        for key in sorted({k for (_t, k) in index}):
            w = index.get(("WT", key))
            m = index.get(("MT", key))
            if not (w and m):
                continue
            pw, pm = w["best"]["pairs"], m["best"]["pairs"]
            fh.write("%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%d\n"
                     % (key, w["best"]["ct"], w["best"]["n_hit"],
                        "" if w["best"]["energy"] is None else "%.2f" % w["best"]["energy"],
                        m["best"]["ct"], m["best"]["n_hit"],
                        "" if m["best"]["energy"] is None else "%.2f" % m["best"]["energy"],
                        ("%.2f" % (m["best"]["energy"] - w["best"]["energy"]))
                        if (w["best"]["energy"] is not None and m["best"]["energy"] is not None) else "",
                        len(pw ^ pm), len(pw & pm)))

    # 4) Log
    with open(os.path.join(out_root, "extract-RNAsnp-sfold-jiaojisubstructure.log.txt"),
              "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(LOG) + "\n")

    log("\n===== Done =====")
    log("Successful tasks: %d (WT %d, MT %d)"
        % (len(results), sum(r["type"] == "WT" for r in results),
           sum(r["type"] == "MT" for r in results)))
    for p in (p1, p2, p3):
        log("  %s" % p)
    return results


def main():
    ap = argparse.ArgumentParser(
        description="Select the single most likely Sfold centroid structure and extract ΔG (RNAsnp probability support)")
    ap.add_argument("--down", default=DEFAULT_DOWN_ROOT, help="sfold-down root directory")
    ap.add_argument("--prob-roots", nargs="*", default=DEFAULT_PROB_ROOTS,
                    help="RNAsnp probability matrix root directories (multiple allowed)")
    ap.add_argument("--out", default=DEFAULT_OUT_ROOT, help="Output root directory")
    ap.add_argument("--type", choices=["WT", "MT", "BOTH"], default="BOTH")
    ap.add_argument("--min-prob", type=float, default=0.0,
                    help="Only count pairings with probability greater than this value (default 0, i.e., p>0 counts as support)")
    ap.add_argument("--no-ensemble", action="store_true",
                    help="Do not include ecentroid.ct as a candidate (default includes it, consistent with old program)")
    ap.add_argument("--dry-run", action="store_true", help="Only count, do not write files")
    args = ap.parse_args()

    types = list(TYPES) if args.type == "BOTH" else [args.type]
    run(args.down, args.prob_roots, args.out, types,
        min_prob=args.min_prob, include_ensemble=not args.no_ensemble,
        dry_run=args.dry_run)


if __name__ == "__main__":
    main()
