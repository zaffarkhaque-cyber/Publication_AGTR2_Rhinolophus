#!/usr/bin/env python3
"""
trace_after_promoter.py - what lies downstream of a promoter whose gene is gone?

FBP2 is lost in hummingbirds, yet ~983 bp of its promoter is still recognisable at 78% identity. This
asks the obvious next question: is there anything left of the coding sequence behind it?

For each species it takes the promoter's own coordinates from PROMOTER.tsv, reads a window of genomic DNA
around them, and searches that window with the human FBP2 protein using tblastn. Hits are chained, then
the chained model is searched back against the complete human RefSeq proteome with blastp. A remnant is
reported as FBP2 only if that search names FBP2 and nothing else at an equal or better score - FBP1 is a
close paralogue and is exactly what a careless search would return instead.

Every in-frame stop codon at least 10 residues inside a hit is reported with its position on the human
protein, its coordinate in the assembly and the codon read back out of the genome.

Apus apus is the positive control: it has FBP2, so an intact coding sequence must be found behind its
promoter. If it is not, the method is wrong and nothing else here means anything.

    python trace_after_promoter.py --reference human.zip --genomes genomes.tsv \\
        --results results --db-dir /path/to/blast_databases \\
        --species Archilochus_colubris,Calypte_anna,Heliangelus_exortis,Apus_apus
"""
import argparse
import csv
import datetime
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import TABLE1 as T                                                       # noqa: E402

GENE = "FBP2"
PARALOG = "FBP1"
WINDOW = 60000           # bp read either side of the promoter hit
STOP_EDGE = 10
EVALUE = "1e-3"


def tblastn_subject(prot_fa, seq, threads):
    """tblastn of a protein against one stretch of DNA, with no database to build."""
    fd, sp = tempfile.mkstemp(suffix=".fna")
    with os.fdopen(fd, "w") as fh:
        fh.write(">window\n%s\n" % seq)
    try:
        out = subprocess.run(
            ["tblastn", "-query", prot_fa, "-subject", sp, "-evalue", EVALUE, "-num_threads", str(threads),
             "-outfmt", "6 qseqid sacc qstart qend sstart send evalue bitscore pident length qlen qseq sseq"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout
    finally:
        os.remove(sp)
    return out


def main():
    ap = argparse.ArgumentParser(description="Trace coding remnants behind a promoter")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--results", default="results", help="directory holding PROMOTER.tsv")
    ap.add_argument("--db-dir", required=True)
    ap.add_argument("--proteome-db-dir", default=None)
    ap.add_argument("--promoter", default="FBP2_promoter_bird",
                    help="which PROMOTER.tsv row supplies the coordinates")
    ap.add_argument("--species", default=None)
    ap.add_argument("--gene", default=GENE)
    ap.add_argument("--window", type=int, default=WINDOW)
    ap.add_argument("--max-models", type=int, default=10,
                    help="chains examined per window; a decayed gene can sit well below its paralogue")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    rdir = T.fix_path(a.results)
    dbdir = T.fix_path(a.db_dir)
    pdbdir = T.fix_path(a.proteome_db_dir) if a.proteome_db_dir else dbdir
    out = T.fix_path(a.out) if a.out else rdir
    start = datetime.datetime.now()

    print("Reading the human reference ...")
    prot, sym = T.load_reference(T.fix_path(a.reference), rdir)
    pdb = T.build_proteome_db(prot, pdbdir)
    acc_q = T.best_protein(a.gene, prot, sym)
    if not acc_q:
        sys.exit("ERROR: %s is not in the reference proteome" % a.gene)
    qfa = os.path.join(out, "trace_%s_query.faa" % a.gene)
    with open(qfa, "w") as fh:
        fh.write(">%s\n%s\n" % (a.gene, prot[acc_q]))
    qlen = len(prot[acc_q])
    print("Query: %s %s, %d aa" % (a.gene, acc_q, qlen))

    # where the promoter was found, taken from the screen's own output
    hits = {}
    with open(os.path.join(rdir, "PROMOTER.tsv"), encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["promoter"] == a.promoter and r["call"] == "DETECTED":
                note = r["anchor"]
                if "best locus" not in note:
                    continue
                loc = note.split("best locus")[1].strip()
                acc, rng = loc.split(":")
                lo, hi = (int(x) for x in rng.split("-"))
                hits[r["species"]] = (acc, lo, hi)
    if not hits:
        sys.exit("ERROR: no DETECTED rows for %s in PROMOTER.tsv" % a.promoter)

    gl = {}
    for ln in open(T.fix_path(a.genomes), encoding="utf-8-sig"):
        if ln.lstrip("﻿").startswith("#") or not ln.strip():
            continue
        f = [x.strip() for x in ln.rstrip("\n").split("\t")]
        if len(f) >= 5:
            gl[f[0]] = dict(acc=f[3], zip=T.fix_path(f[4]))

    want = [s.strip() for s in a.species.split(",")] if a.species else sorted(hits)
    rows, events = [], []
    for sp in want:
        if sp not in hits:
            print("\n==== %s: no %s hit, skipped" % (sp, a.promoter))
            continue
        pacc, plo, phi = hits[sp]
        print("\n==== %s   promoter at %s:%d-%d" % (sp, pacc, plo, phi))
        db, _ = T.build_db(sp, gl[sp]["zip"], dbdir)
        _, seqlen = T.seq_info(db, pacc, {})
        wlo, whi = max(1, plo - a.window), min(seqlen, phi + a.window)
        seq = T.run(["blastdbcmd", "-db", db, "-entry", pacc, "-range", "%d-%d" % (wlo, whi),
                     "-outfmt", "%s"])
        seq = "".join(seq.split()).upper()

        hs = T.parse_hsps(tblastn_subject(qfa, seq, a.threads)).get(a.gene, [])
        if not hs:
            print("  no %s coding sequence anywhere in the %.0f kb window" % (a.gene, (whi - wlo) / 1000.0))
            rows.append([sp, gl[sp]["acc"], "%s:%d-%d" % (pacc, plo, phi), (whi - wlo),
                         "NONE", 0, "-", "-", 0, 0, "-", "-"])
            continue
        # Every chain, not just the best one. FBP1 is an intact neighbour of FBP2 and outscores a decayed
        # FBP2 every time, so taking chains[0] reports the paralogue and declares the gene absent - the
        # same failure this project's other screens were built to avoid.
        chains = T.chain_hsps(hs, 100000)[:a.max_models]
        recs = {}
        for i, ch in enumerate(chains):
            recs["m%d" % i] = "".join(h["sseq"].replace("-", "").replace("*", "X")
                                      for h in sorted(ch["hsps"], key=lambda h: h["qs"]))
        back_all = T.blastp_best(recs, pdb, a.threads, sym)
        named_all, keep = [], None
        for i, ch in enumerate(chains):
            b = back_all.get("m%d" % i)
            nm = b[0] if b else "no blastp match"
            good = bool(b) and b[0] == a.gene and not (b[2] and b[2] != a.gene and b[3] >= b[1])
            named_all.append((nm, int(ch["score"]), good))
            if good and keep is None:
                keep = i
        print("  chains examined: " + ", ".join("%s(%d)%s" % (n, s, "*" if g else "")
                                                for n, s, g in named_all))
        c = chains[keep if keep is not None else 0]
        back = back_all.get("m%d" % (keep if keep is not None else 0))
        named = back[0] if back else "no blastp match"
        ok = keep is not None
        cov = T.coverage(c["hsps"], qlen)
        pid = sum(h["pid"] * h["alen"] for h in c["hsps"]) / sum(h["alen"] for h in c["hsps"])
        stops = []
        for h in c["hsps"]:
            s = h["sseq"]
            for k in range(len(s)):
                if s[k] == "*" and STOP_EDGE <= k < len(s) - STOP_EDGE:
                    stops.append((h, k))
        verdict = ("CODING REMNANT" if ok and stops else
                   "INTACT CODING SEQUENCE" if ok else
                   "NOT %s (best blastp match: %s)" % (a.gene, named))
        glo, ghi = wlo + c["lo"] - 1, wlo + c["hi"] - 1
        print("  %-26s %.0f%% of %s, %.1f%% identity, %d piece(s), %d in-frame stop(s)"
              % (verdict, 100 * cov, a.gene, pid, len(c["hsps"]), len(stops)))
        print("  window %s:%d-%d, model at %s:%d-%d" % (pacc, wlo, whi, pacc, glo, ghi))
        rows.append([sp, gl[sp]["acc"], "%s:%d-%d" % (pacc, plo, phi), (whi - wlo), verdict,
                     len(c["hsps"]), "%.0f%%" % (100 * cov), "%.1f" % pid, len(stops),
                     int(c["score"]), named, "%s:%d-%d" % (pacc, glo, ghi)])
        for h, k in stops:
            qpos = h["qs"] + len(h["qseq"][:k].replace("-", ""))
            events.append([sp, a.gene, qpos, pacc, "see model coordinates", h["sseq"][max(0, k - 5):k + 6]])

    hdr = ["species", "assembly", "promoter locus", "window bp", "verdict", "pieces",
           "coverage of " + a.gene, "identity %", "in-frame stops", "bits",
           "blastp names the model", "model location"]
    p = os.path.join(out, "TRACE_%s.tsv" % a.gene)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(hdr)
        w.writerows(rows)
    p2 = os.path.join(out, "TRACE_%s_stops.tsv" % a.gene)
    with open(p2, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "gene", "position on human protein", "scaffold", "note", "context"])
        w.writerows(events)
    print("\nWrote %s and %s" % (p, p2))
    print("THIS RUN: trace_after_promoter.py md5 %s, finished %s"
          % (T.md5(os.path.abspath(__file__)), datetime.datetime.now().isoformat(timespec="seconds")))
    print("started %s" % start.isoformat(timespec="seconds"))


if __name__ == "__main__":
    main()
