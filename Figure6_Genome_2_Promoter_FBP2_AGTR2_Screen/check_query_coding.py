#!/usr/bin/env python3
"""
check_query_coding.py - does a promoter query contain coding sequence?

A promoter window taken upstream of a gene model can overrun into the gene itself, or into a neighbouring
one. If it does, part of what the screen measures is exonic rather than regulatory, and a "promoter
retained / promoter lost" result is partly a statement about coding sequence.

The test is direct: translate each query in all six frames and search it against the complete human
RefSeq proteome with blastx. A clean promoter returns nothing. Anything that returns a protein is
reported with the gene, the E-value and WHICH BASES of the query are involved, so the query can be
trimmed to the non-coding part if needed.

    python check_query_coding.py --reference human.zip --queries queries.tsv --results results \\
        --db-dir /path/to/blast_databases
"""
import argparse
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import TABLE1 as T                                                       # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Check promoter queries for coding sequence")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--results", default="results")
    ap.add_argument("--db-dir", required=True)
    ap.add_argument("--proteome-db-dir", default=None)
    ap.add_argument("--evalue", default="1e-3")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    rdir = T.fix_path(a.results)
    pdbdir = T.fix_path(a.proteome_db_dir) if a.proteome_db_dir else T.fix_path(a.db_dir)
    print("Reading the human reference ...")
    prot, sym = T.load_reference(T.fix_path(a.reference), rdir)
    pdb = T.build_proteome_db(prot, pdbdir)

    rows = []
    for ln in open(T.fix_path(a.queries), encoding="utf-8-sig"):
        if ln.lstrip("﻿").startswith("#") or not ln.strip():
            continue
        f = [x.strip() for x in ln.rstrip("\n").split("\t")]
        if len(f) < 4:
            continue
        label, src = f[0], f[3]
        path = (os.path.join(HERE, src) if not src.startswith("{") else
                os.path.join(rdir, "queries", "%s.fna" % label))
        path = T.fix_path(path)
        if not os.path.exists(path):
            print("  %-30s query file not found: %s" % (label, path))
            continue
        seqs = T.read_fasta(path)
        qlen = len(list(seqs.values())[0]) if seqs else 0
        out = subprocess.run(
            ["blastx", "-query", path, "-db", pdb, "-evalue", a.evalue, "-num_threads", str(a.threads),
             "-max_target_seqs", "20", "-outfmt", "6 sseqid qstart qend pident length evalue bitscore"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout
        hits = []
        for l in out.splitlines():
            g = l.split("\t")
            if len(g) < 7:
                continue
            hits.append(dict(gene=sym.get(g[0], g[0]), qs=int(g[1]), qe=int(g[2]), pid=float(g[3]),
                             alen=int(g[4]), ev=float(g[5]), bits=float(g[6])))
        if not hits:
            print("  %-30s %4d bp  CLEAN - no protein aligns to this query" % (label, qlen))
            rows.append([label, qlen, "CLEAN", "-", 0, 0, "-", "-"])
            continue
        best = {}
        for h in hits:
            b = best.get(h["gene"])
            if b is None or h["bits"] > b["bits"]:
                best[h["gene"]] = h
        covered = set()
        for h in hits:
            covered.update(range(min(h["qs"], h["qe"]), max(h["qs"], h["qe"]) + 1))
        top = sorted(best.values(), key=lambda h: -h["bits"])[0]
        print("  %-30s %4d bp  CODING: %d bp (%.0f%%) aligns to protein; best %s E=%.1g"
              % (label, qlen, len(covered), 100.0 * len(covered) / max(1, qlen), top["gene"], top["ev"]))
        for h in sorted(best.values(), key=lambda h: -h["bits"])[:4]:
            print("        %-12s query %d-%d, %.0f%% id over %d aa, E=%.1g"
                  % (h["gene"], h["qs"], h["qe"], h["pid"], h["alen"], h["ev"]))
        rows.append([label, qlen, "CODING", top["gene"], len(covered),
                     round(100.0 * len(covered) / max(1, qlen)), "%d-%d" % (top["qs"], top["qe"]),
                     "%.1g" % top["ev"]])

    p = os.path.join(rdir, "QUERY_CODING_CHECK.tsv")
    with open(p, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["query", "query bp", "verdict", "best protein", "coding bp in query",
                    "% of query coding", "query range", "E-value"])
        w.writerows(rows)
    print("\nWrote %s" % p)
    print("A query reported CODING overlaps a gene. Trim it to the non-coding part, or state in the")
    print("methods that the window includes exonic sequence - it changes what the measurement means.")


if __name__ == "__main__":
    main()
