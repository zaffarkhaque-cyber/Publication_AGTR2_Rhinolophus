#!/usr/bin/env python3
"""
INTERGENIC.py - dual-bait scan of the PLS3-SLC6A14 intergenic region for any remaining AGTR2 sequence.

Question: in a genome where AGTR2 is not found, is there a decayed REMNANT (pseudogene) left between the
flanking genes, or is the locus deleted outright?

Method:
  1. The PLS3-SLC6A14 interval is located in each assembly by the TABLE1 method (tblastn + collinear chaining +
     reciprocal-best-hit orthology against the complete human proteome) and its DNA is read straight from the
     assembly. Nothing is downloaded from a browser, so every interval has an assembly accession and coordinates.
  2. Baits are the AGTR2 proteins of the bait species, taken from their own assemblies by the same method (so the
     bait has provenance too, and is not a pasted sequence).
  3. Each bait is searched against each interval with tblastn (protein vs translated DNA). tblastn is used rather
     than blastn because a remnant decaying for tens of millions of years loses nucleotide identity long before it
     loses protein-level similarity - the sensitive direction for a claim of absence.
  4. NOISE FLOOR, measured not assumed: each interval is shuffled --shuffles times and rescanned. The highest score
     obtained from shuffled sequence is that species' empirical floor.
  5. DETECTION LIMIT, so that absence means something: the bait's own coding DNA is mutated at increasing rates,
     embedded in shuffled interval sequence, and put through the same scan. This gives the most diverged remnant
     the scan can still find, i.e. what a negative result excludes.

Controls, printed first and written to INTERGENIC_controls.txt:
  INTACT   species whose AGTR2 is intact must score far above their noise floor (the scan works).
  REMNANT  species whose AGTR2 is a pseudogene must still be detected (the scan finds decayed sequence).
           Without this control a negative result in another species means nothing.
  FLOOR    every species' observed floor from shuffled sequence.
  N        N content of each interval, so an absence cannot be confused with an assembly gap.

Usage:
  python INTERGENIC.py --reference human.zip --genomes genomes.tsv --out results_intergenic
        [--db-dir DIR] [--proteome-db-dir DIR] [--bait-species Myotis_myotis,Apus_apus]
        [--species A,B,C] [--shuffles 20] [--threads 8]
Requires: Python 3.8+, NCBI BLAST+, matplotlib (for the figure only), and TABLE1.py in the same folder.
"""
import argparse
import datetime
import importlib.util
import os
import random
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
# TABLE1.py supplies the gene-finding method (tblastn + chaining + reciprocal-best-hit orthology). It may sit
# beside this script or in a sibling TABLE1 folder, so that each figure can have its own directory.
_CANDIDATES = [os.path.join(HERE, "TABLE1.py"),
               os.path.join(os.path.dirname(HERE), "TABLE1", "TABLE1.py"),
               os.path.join(os.path.dirname(HERE), "TABLE1.py")]
_T1 = next((p for p in _CANDIDATES if os.path.exists(p)), None)
if _T1 is None:
    sys.exit("ERROR: TABLE1.py not found. Looked in:\n  " + "\n  ".join(_CANDIDATES) +
             "\nPut TABLE1.py beside this script or in a sibling TABLE1 folder.")
_spec = importlib.util.spec_from_file_location("TABLE1", _T1)
T = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(T)

BAIT_GENE = "AGTR2"
DIVERGENCE = [0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
EVAL_MAX = 1e-3     # a remnant must also be significant in its own right, not just above the shuffled floor


def detected(m, fl):
    """A remnant is called only if it beats the shuffled-sequence floor AND is significant on its own.
    Either test alone lets through chance alignments in a long interval."""
    return m["nident"] > fl and m["length"] > 0 and m["evalue"] < EVAL_MAX


def dna_of(db, acc, lo, hi, strand="+"):
    s = T.run(["blastdbcmd", "-db", db, "-entry", acc, "-range", "%d-%d" % (lo, hi),
               "-strand", "plus" if strand == "+" else "minus", "-outfmt", "%s"])
    return "".join(s.split()).upper()


BACK = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
        "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT",
        "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}


def backtranslate(prot):
    """Coding DNA for a bait supplied as protein only, for the detection-limit simulation. One fixed codon per
    residue: the simulation measures how much DIVERGENCE the search tolerates, which does not depend on which
    synonymous codon the remnant started from."""
    return "".join(BACK.get(c, "NNN") for c in prot)


def shuffled(seq, rng):
    """Mononucleotide shuffle: same length and base composition, no real sequence."""
    b = list(seq)
    rng.shuffle(b)
    return "".join(b)


def mutate(seq, rate, rng):
    """Random substitutions at the given rate (N positions are left alone)."""
    b = list(seq)
    for i, c in enumerate(b):
        if c in "ACGT" and rng.random() < rate:
            b[i] = rng.choice([x for x in "ACGT" if x != c])
    return "".join(b)


def scan(baits, subjects, threads):
    """tblastn every bait against every subject. subjects: {name: dna}. Returns {(bait, name): metrics}."""
    if not subjects:
        return {}
    fd, sub = tempfile.mkstemp(suffix=".fna")
    with os.fdopen(fd, "w") as fh:
        for n, s in subjects.items():
            fh.write(">%s\n%s\n" % (n, s))
    fd, qry = tempfile.mkstemp(suffix=".faa")
    with os.fdopen(fd, "w") as fh:
        for n, s in baits.items():
            fh.write(">%s\n%s\n" % (n, s))
    try:
        out = T.run(["tblastn", "-query", qry, "-subject", sub, "-evalue", "10", "-outfmt",
                     "6 qseqid sseqid nident bitscore evalue pident length sstart send",
                     "-max_hsps", "1", "-num_threads", str(threads)])
    finally:
        os.remove(sub)
        os.remove(qry)
    best = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 9:
            continue
        k = (f[0], f[1])
        m = dict(nident=int(f[2]), bits=float(f[3]), evalue=float(f[4]), pident=float(f[5]),
                 length=int(f[6]), sstart=int(f[7]), send=int(f[8]))
        if k not in best or m["nident"] > best[k]["nident"]:
            best[k] = m
    return best


def blank():
    return dict(nident=0, bits=0.0, evalue=float("nan"), pident=0.0, length=0, sstart=0, send=0)


# ============================================================================ main
def main():
    ap = argparse.ArgumentParser(description="PLS3-SLC6A14 intergenic scan for AGTR2 remnants")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--db-dir", default=None)
    ap.add_argument("--proteome-db-dir", default=None)
    ap.add_argument("--bait-species", default="Myotis_myotis,Apus_apus",
                    help="species whose AGTR2 proteins are used as baits (must be in the genome list). Ignored for "
                         "the search itself if --baits is given, but still used for the detection-limit simulation.")
    ap.add_argument("--baits", default=None,
                    help="YOUR OWN bait proteins: a FASTA file. Each record becomes one bar group in the figure, "
                         "labelled by its header. Use this to keep baits you downloaded and checked yourself; the "
                         "run log records each one's length and md5 so the figure is tied to the exact sequences.")
    ap.add_argument("--species", default=None, help="comma-separated subset to scan and plot (default: all)")
    ap.add_argument("--shuffles", type=int, default=20, help="shuffled replicates per interval for the noise floor")
    ap.add_argument("--seed", type=int, default=20261117)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--blast-bin", default=None)
    ap.add_argument("--max-intron", type=int, default=1000000)
    ap.add_argument("--min-coverage", type=float, default=0.60)
    a = ap.parse_args()
    start = datetime.datetime.now()
    if a.blast_bin:
        os.environ["PATH"] = T.fix_path(a.blast_bin) + os.pathsep + os.environ["PATH"]
    for t in ("tblastn", "blastp", "makeblastdb", "blastdbcmd"):
        if not __import__("shutil").which(t):
            sys.exit("ERROR: %s not found (install NCBI BLAST+ or give --blast-bin)" % t)
    # check every input file up front: locating the intervals takes minutes, so fail before that, not after
    for lbl, p in (("--reference", a.reference), ("--genomes", a.genomes), ("--baits", a.baits)):
        if p and not os.path.exists(T.fix_path(p)):
            sys.exit("ERROR: %s file not found: %s%s" % (lbl, p, "\n(omit --baits to take the baits from the "
                     "assemblies of --bait-species instead)" if lbl == "--baits" else ""))
    out = T.fix_path(a.out)
    dbdir = T.fix_path(a.db_dir) if a.db_dir else os.path.join(out, "blast_dbs")
    for d in (out, dbdir, os.path.join(out, "raw_tblastn")):
        os.makedirs(d, exist_ok=True)
    ref = T.fix_path(a.reference)
    rng = random.Random(a.seed)

    # ---------------------------------------------------------------- reference and queries (TABLE1 method)
    prot, sym = T.load_reference(ref, out)
    queries, qacc = {}, {}
    for g in T.GENES:
        p = T.best_protein(g, prot, sym)
        if not p:
            sys.exit("ERROR: %s not in the reference annotation" % g)
        queries[g], qacc[g] = prot[p], p
    qfa = os.path.join(out, "INTERGENIC_queries.faa")
    with open(qfa, "w") as fh:
        for g in T.GENES:
            fh.write(">%s\n%s\n" % (g, queries[g]))
    qtag = T.md5(qfa)[:10]
    pdbdir = os.path.expanduser(T.fix_path(a.proteome_db_dir)) if a.proteome_db_dir else dbdir
    os.makedirs(pdbdir, exist_ok=True)
    pdb = T.build_proteome_db(prot, pdbdir)

    genomes = []
    for ln in open(T.fix_path(a.genomes), encoding="utf-8"):
        if ln.strip() and not ln.startswith("#"):
            f = [x.strip() for x in ln.rstrip("\n").split("\t")]
            acc, zp = "", f[3]
            if len(f) >= 5:
                acc, zp = (f[3], f[4]) if re.match(r"GC[AF]_\d+\.\d+$", f[3]) else (f[4], f[3])
            genomes.append(dict(sp=f[0], group=f[1], expect=f[2], acc=acc, zip=T.fix_path(zp)))
    wantlist = [x.strip() for x in a.species.split(",")] if a.species else None   # order is kept: it is the figure's
    want = set(wantlist) if wantlist else None
    baits_wanted = [x.strip() for x in a.bait_species.split(",")]
    todo = [g for g in genomes if want is None or g["sp"] in want or g["sp"] in baits_wanted]
    unknown = [s for s in (wantlist or []) if s not in {g["sp"] for g in genomes}]
    if unknown:
        sys.exit("ERROR: --species names not in %s: %s" % (os.path.basename(a.genomes), ", ".join(unknown)))

    # ---------------------------------------------------------------- locate the interval in each assembly
    print("Locating PLS3-SLC6A14 intervals (TABLE1 method)\n")
    iv, baitseq, baitdna, prov, status = {}, {}, {}, {}, {}
    for g in todo:
        sp = g["sp"]
        if not os.path.exists(g["zip"]):
            status[sp] = "SKIPPED - file not found"
            print("  %-28s %s" % (sp, status[sp]))
            continue
        try:
            db, prov[sp] = T.build_db(sp, g["zip"], dbdir)
        except Exception as e:
            status[sp] = "SKIPPED - %s" % str(e).splitlines()[0]
            print("  %-28s %s" % (sp, status[sp]))
            continue
        r = T.evaluate(sp, db, queries, qfa, qtag, pdb, sym, out, a)
        t2 = r.get("agtr2")
        if sp in baits_wanted:
            if not t2:
                sys.exit("ERROR: bait species %s has no AGTR2 model to take a bait from" % sp)
            baitseq[sp] = "".join(h["sseq"].replace("-", "").replace("*", "X")
                                  for h in sorted(t2["hsps"], key=lambda h: h["qs"]))
            baitdna[sp] = dna_of(db, t2["acc"], t2["lo"], t2["hi"], t2["strand"])
            print("  %-28s BAIT SOURCE  AGTR2 %d aa, locus %d bp" % (sp, len(baitseq[sp]), len(baitdna[sp])))
        info = r.get("interval")
        if not info or "len" not in info:
            status[sp] = ("no interval - PLS3 and SLC6A14 are not on one scaffold" if not info
                          else "interval " + info.get("note", "not read"))
            print("  %-28s %s" % (sp, status[sp]))
            continue
        # where the genome-wide AGTR2 model sits, and how far it is from this interval, so that a zero bar can
        # never be confused with absence from the genome
        where, gap = "-", "-"
        if t2:
            where = "%s:%d-%d" % (t2["acc"], t2["lo"], t2["hi"])
            if t2["acc"] == info["acc"]:
                gap = 0 if info["lo"] <= t2["lo"] and t2["hi"] <= info["hi"] else \
                    (info["lo"] - t2["hi"] if t2["hi"] < info["lo"] else t2["lo"] - info["hi"])
            else:
                gap = "different scaffold"
        iv[sp] = dict(db=db, acc=info["acc"], lo=info["lo"], hi=info["hi"], bp=info["len"],
                      n=info["n"], nrun=info["nrun"], inside=info["agtr2_inside"],
                      call=r["call"], where=where, gap=gap,
                      dna=dna_of(db, info["acc"], info["lo"], info["hi"]))
        status[sp] = "ok"
        print("  %-28s %s:%d-%d  %s bp, %d N (longest %d), AGTR2 call %s" % (
            sp, info["acc"], info["lo"], info["hi"], format(info["len"], ","), info["n"], info["nrun"], r["call"]))
    missing = [s for s in baits_wanted if s not in baitseq]
    if missing:
        sys.exit("ERROR: bait species not available: %s" % ", ".join(missing))

    # YOUR OWN baits replace the extracted ones for the search; the extracted locus DNA is kept for the
    # detection-limit simulation, which needs coding DNA rather than protein.
    simdna = dict(baitdna)
    baitinfo = {b: "extracted from the %s assembly by the TABLE1 method" % b for b in baitseq}
    if a.baits:
        bf = T.fix_path(a.baits)
        supplied, cur = {}, None
        for ln in open(bf, encoding="utf-8", errors="replace"):
            ln = ln.rstrip("\n")
            if ln.startswith(">"):
                cur = ln[1:].strip()
                supplied[cur] = ""
            elif cur and ln.strip():
                supplied[cur] += "".join(ln.split()).upper()
        if not supplied:
            sys.exit("ERROR: no sequences found in %s" % bf)
        baitseq = supplied
        baitinfo = {b: "supplied by the user in %s (%d aa, md5 %s)"
                    % (os.path.basename(bf), len(s), T.md5_text(s) if hasattr(T, "md5_text")
                       else __import__("hashlib").md5(s.encode()).hexdigest()) for b, s in supplied.items()}
        print("\nBaits taken from %s:" % bf)
        for b, s in supplied.items():
            print("  %-44s %4d aa" % (b, len(s)))
        print("  (detection limit simulated on the %s AGTR2 locus DNA)" % ", ".join(simdna))

    # ---------------------------------------------------------------- scan, floor, detection limit
    print("\nScanning (tblastn, %d shuffled replicates per interval)" % a.shuffles)
    real = scan(baitseq, {s: iv[s]["dna"] for s in iv}, a.threads)
    floor, fstat = {}, {}
    for s in iv:
        subs = {"shuf_%d" % i: shuffled(iv[s]["dna"], rng) for i in range(a.shuffles)}
        sh = scan(baitseq, subs, a.threads)
        floor[s], fstat[s] = {}, {}
        for b in baitseq:
            v = [sh[(b, k)]["nident"] if (b, k) in sh else 0 for k in subs]
            mu = sum(v) / len(v)
            sd = (sum((x - mu) ** 2 for x in v) / max(1, len(v) - 1)) ** 0.5
            # the observed maximum of a few shuffles underestimates the tail, so take the higher of the
            # observed maximum and a 3-SD upper bound. A score must beat THAT to count as a remnant.
            floor[s][b] = int(round(max(max(v), mu + 3 * sd)))
            fstat[s][b] = (mu, sd, max(v))
        print("  %-28s floor %s" % (s, ", ".join("%s %d (shuffled mean %.1f, sd %.1f, max %d)"
              % (b.split("_")[0], floor[s][b], *fstat[s][b]) for b in baitseq)))

    print("\nDetection limit (bait coding DNA mutated and embedded in shuffled interval sequence)")
    det = {}
    for b in baitseq:
        src = simdna.get(b) or backtranslate(baitseq[b])     # own locus DNA if we have it, else from the protein
        bg = shuffled(iv[list(iv)[0]]["dna"], rng)[:200000]
        subs = {}
        for d in DIVERGENCE:
            m = mutate(src, d, rng)
            half = len(bg) // 2
            subs["div_%02d" % int(d * 100)] = bg[:half] + m + bg[half:]
        sc = scan({b: baitseq[b]}, subs, a.threads)
        gl = max(floor[s][b] for s in floor) if floor else 0
        det[b] = []
        for d in DIVERGENCE:
            m = sc.get((b, "div_%02d" % int(d * 100)), blank())
            det[b].append((d, m["nident"], m["bits"], m["evalue"], detected(m, gl)))
        ok = [d for d, n, bi, e, p in det[b] if p]
        print("  %-28s detected up to %.0f%% divergence (global floor %d identical residues)"
              % (b, 100 * max(ok) if ok else 0, gl))

    # ---------------------------------------------------------------- outputs
    # the figure follows the order given on the command line; without --species, the genome-list order
    order = [s for s in wantlist if s in iv] if wantlist else [g["sp"] for g in todo if g["sp"] in iv]
    baits = list(baitseq)
    with open(os.path.join(out, "INTERGENIC.tsv"), "w") as fh:
        fh.write("species\tassembly accession\tscaffold\tinterval start\tinterval end\t"
                 "interval bp\tN in interval\tlongest N run\tAGTR2 in this interval\t"
                 + "\t".join("%s: identical residues\t%s: bit score\t%s: E-value\t%s: %% identity\t"
                             "%s: aligned aa\t%s: noise floor\t%s: shuffled mean\t%s: shuffled SD\t"
                             "%s: shuffled max\t%s: above floor" % ((b,) * 10) for b in baits) + "\n")
        for s in order:
            v = iv[s]
            row = [s, T.accession(prov.get(s, "")), v["acc"], v["lo"], v["hi"], v["bp"], v["n"], v["nrun"],
                   v["inside"]]
            for b in baits:
                m = real.get((b, s), blank())
                mu, sd, mx = fstat[s][b]
                row += [m["nident"], "%.0f" % m["bits"], "%.1g" % m["evalue"] if m["length"] else "-",
                        "%.1f" % m["pident"], m["length"], floor[s][b], "%.1f" % mu, "%.1f" % sd, mx,
                        "yes" if detected(m, floor[s][b]) else "no"]
            fh.write("\t".join(map(str, row)) + "\n")
        for g in todo:
            if g["sp"] not in iv and g["sp"] in status:
                fh.write("%s\t%s\t-\t-\t-\t-\t-\t-\t-\t%s%s\n" % (
                    g["sp"], T.accession(prov.get(g["sp"], "")), status[g["sp"]], "\t-" * (10 * len(baits) - 1)))

    with open(os.path.join(out, "INTERGENIC_detection_limit.tsv"), "w") as fh:
        fh.write("bait\tdivergence from bait (%)\tidentical residues\tbit score\tE-value\tabove global noise floor\n")
        for b in baits:
            for d, n, bi, e, p in det[b]:
                fh.write("%s\t%.0f\t%d\t%.0f\t%s\t%s\n" % (b, 100 * d, n, bi, "%.1g" % e if n else "-",
                                                           "yes" if p else "no"))

    # ---------------------------------------------------------------- controls
    L = ["CONTROLS"]
    exp = {g["sp"]: g["expect"] for g in genomes}
    # Everything below is about THIS interval only. A species counts as a control only if its AGTR2 lies
    # between the flanking genes; every other species is simply one whose interval was scanned.
    inside = lambda s: str(iv[s]["inside"]).startswith("yes")
    intact = [s for s in order if iv[s]["call"] == "FOUND" and inside(s)]
    remnant = [s for s in order if iv[s]["call"] == "FOUND_DISRUPTED" and inside(s)]
    gone = [s for s in order if s not in intact and s not in remnant]
    for name, group, need in (("INTACT  (intact AGTR2 lies in this interval; the scan must find it)", intact, True),
                              ("REMNANT (an AGTR2 pseudogene lies in this interval; the scan must still find it)",
                               remnant, True)):
        L.append("  " + name)
        if not group:
            L.append("      none in this run" + ("   <-- WITHOUT THIS CONTROL A NEGATIVE MEANS NOTHING"
                                                 if "REMNANT" in name else ""))
        for s in group:
            hit = {b: detected(real.get((b, s), blank()), floor[s][b]) for b in baits}
            L.append("      %-28s %s  %s" % (s, ", ".join("%s %d (floor %d)" % (b.split("_")[0],
                     real.get((b, s), blank())["nident"], floor[s][b]) for b in baits),
                     "PASS" if any(hit.values()) == need else "FAIL"))
    L.append("  SCANNED, NO AGTR2 IN THE INTERVAL")
    for s in gone:
        above = [b for b in baits if detected(real.get((b, s), blank()), floor[s][b])]
        L.append("      %-28s %s  %s" % (s, ", ".join("%s %d (floor %d)" % (b.split("_")[0],
                 real.get((b, s), blank())["nident"], floor[s][b]) for b in baits),
                 "no remnant above floor" if not above else "REMNANT DETECTED by " + ", ".join(above)))
    L.append("  DETECTION LIMIT (what a negative excludes)")
    for b in baits:
        ok = [d for d, n, bi, e, p in det[b] if p]
        L.append("      %-28s a remnant up to %.0f%% diverged from this bait would have been detected"
                 % (b, 100 * max(ok) if ok else 0))
    L.append("  ASSEMBLY GAPS (an absence inside a gap would not count)")
    for s in order:
        L.append("      %-28s %s bp interval, %d N, longest run %d" % (s, format(iv[s]["bp"], ","), iv[s]["n"],
                                                                       iv[s]["nrun"]))
    txt = "\n".join(L)
    print("\n" + txt)
    open(os.path.join(out, "INTERGENIC_controls.txt"), "w").write(txt + "\n")

    # ---------------------------------------------------------------- methods + log
    tver = T.run(["tblastn", "-version"]).splitlines()[0]
    gl = {b: max(floor[s][b] for s in floor) for b in baits}
    lim = {b: 100 * max([d for d, n, bi, e, p in det[b] if p] or [0]) for b in baits}
    meth = (
        "METHODS (intergenic scan)\nGenerated by INTERGENIC.py from this run.\n\n"
        "The PLS3-SLC6A14 intergenic region was located in each assembly without using its annotation: the human "
        "RefSeq proteins PLS3 (%s) and SLC6A14 (%s) were searched with tblastn (%s), hits were chained into "
        "collinear gene models, and orthology was assigned by reciprocal best hit against the complete human "
        "proteome. Where both flanking orthologues lay on one scaffold, the DNA between them was read directly "
        "from the assembly, and its length, total N content and longest run of N recorded.\n\n"
        "Baits.\n%s\nEach bait was searched against each interval with tblastn (protein against translated DNA), "
        "which retains sensitivity to diverged remnants that nucleotide comparison loses. The score plotted is the "
        "number of identical amino acids in the best alignment; bit scores, E-values and alignment lengths are in "
        "INTERGENIC.tsv.\n\n"
        "Criterion for a remnant. A hit was counted as a remnant only if it both exceeded the noise floor defined "
        "below AND had E < %g; either test on its own admits chance alignments in a long interval. The noise floor "
        "was measured rather than assumed: each interval was shuffled %d times, preserving its length and base "
        "composition, and rescanned. The floor is the greater of the highest score obtained from shuffled sequence "
        "and the mean plus three standard deviations of those scores, so that a limited number of replicates does "
        "not understate the tail (%s). Per-interval shuffled means, standard deviations and maxima are in "
        "INTERGENIC.tsv.\n\n"
        "Detection limit. Each bait's coding sequence was mutated at rates of %s, embedded in shuffled sequence "
        "and put through the same scan and the same criterion: %s. %s\n\n"
        "Controls. Species whose AGTR2 is intact and species whose AGTR2 is a pseudogene were both required to "
        "meet the criterion; the pseudogene species establish that the scan detects decayed sequence, without "
        "which a negative result would be uninterpretable.\n" % (
            qacc["PLS3"], qacc["SLC6A14"], tver,
            "\n".join("  %s: %s" % (b, baitinfo.get(b, "source not recorded")) for b in baits),
            EVAL_MAX, a.shuffles,
            ", ".join("%s %d identical residues" % (b, gl[b]) for b in baits),
            ", ".join("%.0f%%" % (100 * d) for d in DIVERGENCE[1:]),
            ("a remnant up to %.0f%% diverged from either bait would have been detected" % lim[baits[0]]
             if len(set(lim.values())) == 1 else
             "the most diverged remnant still detected was " + ", ".join("%s %.0f%%" % (b, lim[b]) for b in baits)),
            "Coding sequence for this simulation was taken from the AGTR2 locus of %s." % ", ".join(simdna)
            if not a.baits else
            "Because the baits were supplied as protein, the simulated remnant was reverse-translated from each "
            "bait with one fixed codon per residue; the simulation measures tolerated divergence, which does not "
            "depend on the starting synonymous codons."))
    open(os.path.join(out, "INTERGENIC_methods.txt"), "w").write(meth)
    with open(os.path.join(out, "INTERGENIC_run_log.txt"), "w") as fh:
        fh.write("INTERGENIC.py md5 %s\nTABLE1.py md5 %s (%s)\nstarted %s finished %s\n%s\n" % (
            T.md5(os.path.abspath(__file__)), T.md5(_T1), _T1,
            start.isoformat(timespec="seconds"), datetime.datetime.now().isoformat(timespec="seconds"), tver))
        fh.write("command: %s\nseed %d, shuffles %d\n\n" % (" ".join(sys.argv), a.seed, a.shuffles))
        for s in order:
            fh.write("%-28s %s %s:%d-%d\n" % (s, T.accession(prov.get(s, "")), iv[s]["acc"], iv[s]["lo"], iv[s]["hi"]))
        fh.write("\n" + txt + "\n")

    # ---------------------------------------------------------------- figure
    try:
        figure(out, order, baits, iv, real, floor, [g for g in genomes])
    except ImportError:
        print("\n(matplotlib not installed - figure skipped; INTERGENIC.tsv has the numbers)")
    print("\nResults in %s: INTERGENIC.tsv, INTERGENIC_controls.txt, INTERGENIC_methods.txt, "
          "INTERGENIC_detection_limit.tsv, INTERGENIC_run_log.txt, INTERGENIC_scan.png/.pdf" % out)


DARKRED = "#8B1A1A"
BARCOL = ["#2E7DF7", "#FFD21E", "#21918c", "#e6550d"]


def figure(out, order, baits, iv, real, floor, genomes):
    """The published figure: grouped bars, one group per species, one bar per bait, value labels, a dashed
    noise-floor line, and the Chiroptera / Aves split. Species appear in the order given on the command line."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cls = {g["sp"]: ("Aves" if g["group"] in ("Hummingbirds", "Other birds") else "Chiroptera") for g in genomes}
    n, nb = len(order), max(1, len(baits))
    w = 0.74 / nb
    fig, ax = plt.subplots(figsize=(max(9.5, 1.45 * n + 2.5), 6.6))

    allys = []
    for j, b in enumerate(baits):
        xs = [i + (j - (nb - 1) / 2) * w for i in range(n)]
        ys = [real.get((b, s), blank())["nident"] for s in order]
        allys += ys
        ax.bar(xs, ys, width=w, color=BARCOL[j % len(BARCOL)], edgecolor="black", linewidth=0.8,
               label="$\\it{%s}$ bait" % b.replace("_", "\\ "), zorder=3)
    top = max(allys + [10])
    ymax = 50 * (int(top * 1.12 / 50) + 1)
    for j, b in enumerate(baits):
        xs = [i + (j - (nb - 1) / 2) * w for i in range(n)]
        for x, y in zip(xs, [real.get((b, s), blank())["nident"] for s in order]):
            ax.text(x, y + ymax * 0.012, str(y), ha="center", va="bottom", fontsize=8.5, color="#333333", zorder=4)

    gl = max(max(floor[s].values()) for s in order) if order else 0
    ax.axhline(gl, ls="--", color="#C00000", lw=1.4, zorder=2)
    ax.text(n - 0.45, gl + ymax * 0.012, "Noise floor", color="#C00000", fontsize=9, style="italic",
            va="bottom", ha="right", zorder=4)

    first_bird = next((i for i, s in enumerate(order) if cls.get(s) == "Aves"), None)
    if first_bird:
        ax.axvline(first_bird - 0.5, color="#BBBBBB", ls=":", lw=1.3, zorder=1)
        ax.text((first_bird - 1) / 2, -0.165, "Chiroptera", transform=ax.get_xaxis_transform(),
                ha="center", fontsize=10.5, fontweight="bold", color=DARKRED)
        ax.text((first_bird + n - 1) / 2, -0.165, "Aves", transform=ax.get_xaxis_transform(),
                ha="center", fontsize=10.5, fontweight="bold", color=DARKRED)

    ax.set_xticks(range(n))
    ax.set_xticklabels(["%s. %s" % (s.split("_")[0][0], s.split("_", 1)[1].replace("_", " ")) for s in order],
                       fontsize=10.5, fontweight="bold", fontstyle="italic")
    ax.set_ylim(0, ymax)
    ax.set_xlim(-0.65, n - 0.35)
    ax.set_ylabel("Homologous Identity Score\n(identical residues in ungapped window)",
                  fontsize=11, fontweight="bold", color=DARKRED)
    ax.set_title("AGTR2 Intergenic Dual-Bait Scan", fontsize=14, fontweight="bold", color=DARKRED, pad=16)
    ax.grid(axis="y", ls="--", alpha=0.45, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "INTERGENIC_scan.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(out, "INTERGENIC_scan.pdf"), bbox_inches="tight")


if __name__ == "__main__":
    main()
