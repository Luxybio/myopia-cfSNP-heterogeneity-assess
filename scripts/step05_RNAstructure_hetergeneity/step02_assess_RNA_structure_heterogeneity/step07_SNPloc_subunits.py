#!/usr/bin/env python3
# Replace <LOCAL_DATA_ROOT>/<SERVER_DATA_ROOT> with your own paths.
# -*- coding: utf-8 -*-
"""Locate the SNP position within WT/MT substructure units."""

import os
import sys

SEQ_ROOT = r"<LOCAL_DATA_ROOT>\Sfold-seq-build"
LONG_ROOT = r"<LOCAL_DATA_ROOT>\substructure-compute\summary"
OUT_DIR = r"<LOCAL_DATA_ROOT>\Figures-reorganized\FIgures4_revisualise"
TYPES = ("WT", "MT")

LOG = []


def log(msg):
    LOG.append(msg)
    print(msg, flush=True)


def read_seq(path):
    """Read a single FASTA sequence, return the sequence string."""
    seq = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    return "".join(seq)


def find_snp_pos(wt, mt):
    """Compare position by position and return the 1-based differing position; return None if identical."""
    if len(wt) != len(mt):
        return None
    diffs = [i + 1 for i, (a, b) in enumerate(zip(wt, mt)) if a != b]
    return diffs[0] if len(diffs) == 1 else (None if not diffs else diffs)


def load_long(typ):
    """Read {typ}-substructure-long.txt -> {ID: [row dict, ...]}"""
    path = os.path.join(LONG_ROOT, "%s-substructure-long.txt" % typ)
    table = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        idx = {name: k for k, name in enumerate(header)}
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            f = line.split("\t")
            rec = {
                "ID": f[idx["ID"]],
                "class": f[idx["class"]],
                "index": f[idx["index"]],
                "positions": [int(x) for x in f[idx["positions"]].split(",") if x],
            }
            table.setdefault(rec["ID"], []).append(rec)
    return table


def locate_snp_subunit(snp_pos, rows):
    """Find the subunit containing the SNP position in the given structure's subunit list; return the row or None."""
    for r in rows:
        if snp_pos in r["positions"]:
            return r
    return None


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    longs = {t: load_long(t) for t in TYPES}
    log("[INFO] Substructure details loaded: WT %d structures, MT %d structures"
        % (len(longs["WT"]), len(longs["MT"])))

    # jobkey = local sequence file name (without .txt); each exists in both WT and MT directories
    wt_files = sorted(f[:-4] for f in os.listdir(os.path.join(SEQ_ROOT, "WT"))
                      if f.endswith(".txt"))
    log("[INFO] WT local sequences: %d" % len(wt_files))

    results = []
    for jobkey in wt_files:
        wt_seq = read_seq(os.path.join(SEQ_ROOT, "WT", jobkey + ".txt"))
        mt_seq = read_seq(os.path.join(SEQ_ROOT, "MT", jobkey + ".txt"))
        snp_pos = find_snp_pos(wt_seq, mt_seq)

        rec = {"jobkey": jobkey, "snp_pos": snp_pos,
               "wt_seq": wt_seq, "mt_seq": mt_seq,
               "WT": None, "MT": None, "note": ""}

        if snp_pos is None:
            rec["note"] = "WT/MT sequences identical; SNP is outside the window (or cannot be located), no assignment"
            log("[WARN] %s: %s" % (jobkey, rec["note"]))
            results.append(rec)
            continue
        if isinstance(snp_pos, list):
            rec["note"] = "Multiple differences: %s" % snp_pos
            log("[WARN] %s: %s" % (jobkey, rec["note"]))
            results.append(rec)
            continue

        for typ in TYPES:
            # Long-table ID looks like WT_<jobkey>_<centroid>; prefix uniquely matches jobkey
            cand = [iid for iid in longs[typ]
                    if iid.startswith("%s_%s_c" % (typ, jobkey))]
            if not cand:
                log("[WARN] %s %s: corresponding structure ID not found in long table" % (typ, jobkey))
                continue
            rows = longs[typ][cand[0]]
            hit = locate_snp_subunit(snp_pos, rows)
            if hit:
                rec[typ] = {"struct_id": cand[0], "class": hit["class"],
                            "index": hit["index"],
                            "positions": ",".join(map(str, hit["positions"])),
                            "length": len(hit["positions"])}
            else:
                log("[WARN] %s %s: SNP position %s did not fall into any subunit?" % (typ, jobkey, snp_pos))
        results.append(rec)
        w = rec["WT"] or {}
        m = rec["MT"] or {}
        log("[OK] %s SNP@%s  WT:%s(%s)  MT:%s(%s)"
            % (jobkey, snp_pos,
               w.get("class", "NA"), w.get("positions", "NA"),
               m.get("class", "NA"), m.get("positions", "NA")))

    # ---- Output 1: long table ----
    p1 = os.path.join(OUT_DIR, "snp-subunit-location-long.txt")
    with open(p1, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("type\tjobkey\tsnp_pos_in_window\tstruct_id\tsubunit_class\t"
                 "subunit_index\tsubunit_positions\tsubunit_length\tnote\n")
        for r in results:
            for typ in TYPES:
                hit = r[typ]
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    typ, r["jobkey"],
                    "" if r["snp_pos"] is None else r["snp_pos"],
                    hit["struct_id"] if hit else "",
                    hit["class"] if hit else "NA",
                    hit["index"] if hit else "",
                    hit["positions"] if hit else "NA",
                    hit["length"] if hit else "NA",
                    r["note"]))
            if r["WT"] is None and r["MT"] is None and not r["note"]:
                pass

    # ---- Output 2: wide table (WT/MT side by side) ----
    p2 = os.path.join(OUT_DIR, "snp-subunit-location-wide.txt")
    with open(p2, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("jobkey\tsnp_pos_in_window\t"
                 "WT_subunit\tWT_positions\tWT_len\t"
                 "MT_subunit\tMT_positions\tMT_len\t"
                 "same_class\tnote\n")
        for r in results:
            w, m = r["WT"] or {}, r["MT"] or {}
            same = ""
            if w and m:
                same = "1" if w["class"] == m["class"] else "0"
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                r["jobkey"],
                "" if r["snp_pos"] is None else r["snp_pos"],
                w.get("class", "NA"), w.get("positions", "NA"), w.get("length", "NA"),
                m.get("class", "NA"), m.get("positions", "NA"), m.get("length", "NA"),
                same, r["note"]))

    with open(os.path.join(OUT_DIR, "snp-subunit-location.log.txt"), "w",
              encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(LOG) + "\n")

    # ---- Summary ----
    valid = [r for r in results if r["WT"] and r["MT"]]
    diff_class = [r for r in valid if r["WT"]["class"] != r["MT"]["class"]]
    log("\n===== Done =====")
    log("Total jobkeys: %d, SNP position locatable and WT/MT both assigned: %d" % (len(results), len(valid)))
    log("Subunit type containing the SNP changed (WT/MT class differ): %d" % len(diff_class))
    log("  %s" % ", ".join("%s@%s %s->%s" % (r["jobkey"].split("_")[1], r["snp_pos"],
        r["WT"]["class"], r["MT"]["class"]) for r in diff_class))
    log("Output: %s" % p1)
    log("      %s" % p2)
    return results


if __name__ == "__main__":
    main()
