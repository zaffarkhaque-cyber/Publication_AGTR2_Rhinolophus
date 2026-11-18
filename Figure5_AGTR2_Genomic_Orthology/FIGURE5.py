#!/usr/bin/env python3
"""
FIGURE5.py - whole-genome AGTR2 orthology screen.

Question: searching the ENTIRE assembly, not just the AGTR2 locus, is there an AGTR2 orthologue anywhere -
intact, pseudogenised, or translocated - or is the best hit merely a paralogue such as AGTR1?

Method, per genome:
  1. The AGTR2 bait protein is searched against the whole assembly with tblastn.
  2. Hits are CHAINED into collinear gene models (same scaffold and strand, introns <= --max-intron bp).
     Chaining matters here: a pseudogene broken by stop codons and frameshifts is recovered in several pieces,
     and no single piece need outscore one contiguous alignment to a paralogue. Taking the best single HSP
     therefore reports the paralogue and misses the pseudogene - the error this script exists to avoid.
  3. Each model's conceptual translation is searched (blastp) against the COMPLETE human RefSeq proteome.
     The model is AGTR2 only if its best match is AGTR2 itself; otherwise it is reported under the gene it
     actually matches (AGTR1, BDKRB1, CXCR3 ... all of which cross-hit an AGTR2 query at 27-38% identity).
  4. Each species is then called:
       INTACT      an AGTR2 orthologue covering >= --min-coverage of the bait, with no internal stop codon
       PSEUDOGENE  an AGTR2 orthologue carrying internal stop codons and/or recovered in several pieces
       NO_AGTR2    no AGTR2 orthologue anywhere; the best hit is a named paralogue
  5. The PARALOGUE FLOOR is measured, not assumed: for every species the best NON-AGTR2 model is reported
     with the gene it belongs to and its bit score. The dashed line in the figure is the highest of those,
     so the threshold comes from the data rather than from a number chosen by eye.

Controls (printed first, written to FIGURE5_controls.txt):
  PRESENT  species whose AGTR2 is reported in the literature must come out INTACT.
  ABSENT   species with a published loss must come out NO_AGTR2.
  AGTR1    AGTR1 is present in every vertebrate, so it must be found in every assembly. A genome where the
           scan cannot find AGTR1 is a genome where the scan is not working, and its result means nothing.

Usage:
  python FIGURE5.py --reference human.zip --genomes genomes.tsv --out results
        [--baits my_baits.faa] [--bait-species Myotis_myotis] [--species A,B,C]
        [--db-dir DIR] [--proteome-db-dir DIR] [--threads 8]
Requires: Python 3.8+, NCBI BLAST+, matplotlib (figure only), and TABLE1.py beside this script.
"""
import argparse
import datetime
import importlib.util
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_CANDIDATES = [os.path.join(HERE, "TABLE1.py"),
               os.path.join(os.path.dirname(HERE), "TABLE1", "TABLE1.py"),
               os.path.join(os.path.dirname(HERE), "TABLE1.py")]
_T1 = next((p for p in _CANDIDATES if os.path.exists(p)), None)
if _T1 is None:
    sys.exit("ERROR: TABLE1.py not found. Looked in:\n  " + "\n  ".join(_CANDIDATES))
_spec = importlib.util.spec_from_file_location("TABLE1", _T1)
T = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(T)

TARGET = "AGTR2"
PARALOG = "AGTR1"          # present in every vertebrate: used as a positive control on the search itself
EVALUE = "1e-3"
STOP_EDGE = 10             # a stop codon this far inside a hit is counted; nearer the ends it may be a junction


def models_for(sp, db, bait, baitname, pdb, sym, out, a):
    """Every chained gene model the bait finds in this assembly, each labelled with the gene it really matches."""
    cache = os.path.join(out, "raw_tblastn", "%s.%s.tsv" % (sp, a._qtag))
    hs = T.parse_hsps(T.tblastn(a._qfa, db, a.threads, cache))
    chains = T.chain_hsps(hs.get(baitname, []), a.max_intron)
    if chains:
        chains = [c for c in chains if c["score"] >= 0.1 * chains[0]["score"]][:a.max_models]
    recs = {}
    for i, c in enumerate(chains):
        c.update(rank=i + 1, id="M%d" % i,
                 cov=T.coverage(c["hsps"], len(bait)),
                 pid=sum(h["pid"] * h["alen"] for h in c["hsps"]) / sum(h["alen"] for h in c["hsps"]),
                 stops=sum(h["sseq"][STOP_EDGE:-STOP_EDGE].count("*") for h in c["hsps"]),
                 evalue=min(h["evalue"] for h in c["hsps"]))
        p = "".join(h["sseq"].replace("-", "").replace("*", "X")
                    for h in sorted(c["hsps"], key=lambda h: h["qs"]))
        if len(p) >= 25:
            recs[c["id"]] = p
    back = T.blastp_best(recs, pdb, a.threads, sym)
    info = {}
    for c in chains:
        b = back.get(c["id"])
        c["gene"] = b[0] if b else "no blastp match"
        # same rule as TABLE1: the best human match must be AGTR2, and no other gene may tie or beat it
        c["is_target"] = bool(b) and b[0] == TARGET and not (b[2] and b[2] != TARGET and b[3] >= b[1])
        c["chrom"], _ = T.seq_info(db, c["acc"], info)
    return chains


def main():
    ap = argparse.ArgumentParser(description="Whole-genome AGTR2 orthology screen (Figure 5)")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--baits", default=None, help="FASTA with the bait protein; the FIRST record is used")
    ap.add_argument("--bait-species", default="Myotis_myotis",
                    help="if --baits is not given, take the bait from this species' assembly")
    ap.add_argument("--species", default=None, help="comma-separated subset, in figure order (default: all)")
    ap.add_argument("--db-dir", default=None)
    ap.add_argument("--proteome-db-dir", default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--blast-bin", default=None)
    ap.add_argument("--max-intron", type=int, default=1000000)
    ap.add_argument("--min-coverage", type=float, default=0.60)
    ap.add_argument("--max-models", type=int, default=25,
                    help="how many of the highest-scoring gene models to identify per assembly. The score of "
                         "the lowest one examined is reported, so the depth of the search is stated, not "
                         "assumed; it must fall below the measured paralogue floor for a negative to mean "
                         "anything.")
    ap.add_argument("--no-figure", action="store_true",
                    help="write the tables and controls only; draw nothing")
    ap.add_argument("--divergent-below", type=float, default=0.0, metavar="PCT",
                    help="shade intact AGTR2 orthologues below this %% identity in a paler colour, and say "
                         "so in the legend. Presentation only: the call and every number are unchanged, and "
                         "the cutoff is printed on the figure so it is never an unexplained judgement. "
                         "0 (the default) draws every intact orthologue the same.")
    a = ap.parse_args()
    start = datetime.datetime.now()
    import shutil as _sh
    if a.blast_bin:
        os.environ["PATH"] = T.fix_path(a.blast_bin) + os.pathsep + os.environ["PATH"]
    for t in ("tblastn", "blastp", "makeblastdb", "blastdbcmd"):
        if _sh.which(t) is None:
            sys.exit("ERROR: %s not found (install NCBI BLAST+ or give --blast-bin)" % t)
    for lbl, p in (("--reference", a.reference), ("--genomes", a.genomes), ("--baits", a.baits)):
        if p and not os.path.exists(T.fix_path(p)):
            sys.exit("ERROR: %s file not found: %s" % (lbl, p))
    out = T.fix_path(a.out)
    dbdir = T.fix_path(a.db_dir) if a.db_dir else os.path.join(out, "blast_dbs")
    for d in (out, dbdir, os.path.join(out, "raw_tblastn")):
        os.makedirs(d, exist_ok=True)

    prot, sym = T.load_reference(T.fix_path(a.reference), out)
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
    wantlist = [x.strip() for x in a.species.split(",")] if a.species else None
    unknown = [s for s in (wantlist or []) if s not in {g["sp"] for g in genomes}]
    if unknown:
        sys.exit("ERROR: --species names not in %s: %s" % (os.path.basename(a.genomes), ", ".join(unknown)))
    todo = [g for g in genomes if wantlist is None or g["sp"] in wantlist or g["sp"] == a.bait_species]

    # ---------------------------------------------------------------- bait
    baitname, bait, baitinfo = None, None, ""
    if a.baits:
        bf = T.fix_path(a.baits)
        cur = None
        for ln in open(bf, encoding="utf-8", errors="replace"):
            ln = ln.rstrip("\n")
            if ln.startswith(">"):
                if bait:
                    break                      # only the first record is used
                baitname, bait = ln[1:].strip(), ""
            elif baitname is not None and ln.strip():
                bait += "".join(ln.split()).upper()
        if not bait:
            sys.exit("ERROR: no sequence found in %s" % bf)
        import hashlib
        baitinfo = "supplied by the user in %s (%d aa, md5 %s)" % (
            os.path.basename(bf), len(bait), hashlib.md5(bait.encode()).hexdigest())
    print("Whole-genome AGTR2 orthology screen\n")

    res, prov, status = {}, {}, {}
    pending = []
    for g in todo:
        if not os.path.exists(g["zip"]):
            status[g["sp"]] = "SKIPPED - file not found: %s" % g["zip"]
            print("  %-28s %s" % (g["sp"], status[g["sp"]]))
            continue
        try:
            db, prov[g["sp"]] = T.build_db(g["sp"], g["zip"], dbdir)
        except Exception as e:
            status[g["sp"]] = "SKIPPED - %s" % str(e).splitlines()[0]
            print("  %-28s %s" % (g["sp"], status[g["sp"]]))
            continue
        pending.append((g, db))

    if bait is None:                               # take the bait from an assembly, by the TABLE1 method
        src = next((x for x in pending if x[0]["sp"] == a.bait_species), None)
        if src is None:
            sys.exit("ERROR: bait species %s is not available; give --baits instead" % a.bait_species)
        q = {gg: prot[T.best_protein(gg, prot, sym)] for gg in T.GENES}
        qf = os.path.join(out, "_tmp_queries.faa")
        with open(qf, "w") as fh:
            for gg in T.GENES:
                fh.write(">%s\n%s\n" % (gg, q[gg]))
        a._qfa, a._qtag = qf, T.md5(qf)[:10]
        r = T.evaluate(a.bait_species, src[1], q, qf, a._qtag, pdb, sym, out, a)
        if not r.get("agtr2"):
            sys.exit("ERROR: no AGTR2 model in %s to take a bait from" % a.bait_species)
        bait = "".join(h["sseq"].replace("-", "").replace("*", "X")
                       for h in sorted(r["agtr2"]["hsps"], key=lambda h: h["qs"]))
        baitname = a.bait_species
        baitinfo = "extracted from the %s assembly (%s) by the TABLE1 method (%d aa)" % (
            a.bait_species, T.accession(prov.get(a.bait_species, "")), len(bait))
    print("Bait: %s - %s\n" % (baitname, baitinfo))

    bfa = os.path.join(out, "FIGURE5_bait.faa")
    with open(bfa, "w") as fh:
        fh.write(">%s\n%s\n" % (baitname, bait))
    a._qfa, a._qtag = bfa, T.md5(bfa)[:10]

    for g, db in pending:
        sp = g["sp"]
        ms = models_for(sp, db, bait, baitname, pdb, sym, out, a)
        tgt = [c for c in ms if c["is_target"]]
        par = [c for c in ms if not c["is_target"]]
        best_t = max(tgt, key=lambda c: c["score"]) if tgt else None
        best_p = max(par, key=lambda c: c["score"]) if par else None
        agtr1 = [c for c in par if c["gene"] == PARALOG]
        if best_t and best_t["cov"] >= a.min_coverage and best_t["stops"] == 0 and len(best_t["hsps"]) == 1:
            call = "INTACT"
        elif best_t:
            call = "PSEUDOGENE"
        else:
            call = "NO_AGTR2"
        res[sp] = dict(call=call, target=best_t, paralog=best_p, models=ms,
                       nmodels=len(ms), lowest=min([c["score"] for c in ms] or [0]),
                       agtr1=max(agtr1, key=lambda c: c["score"]) if agtr1 else None)
        status[sp] = "ok"
        b = best_t or best_p
        print("  %-28s %-11s best hit %-9s %4.0f bits, %4.1f%% id, %3.0f%% cov, %d piece(s), %d stop(s)" % (
            sp, call, (b["gene"] if b else "none"), b["score"] if b else 0, b["pid"] if b else 0,
            100 * b["cov"] if b else 0, len(b["hsps"]) if b else 0, b["stops"] if b else 0))

    order = [s for s in (wantlist or [g["sp"] for g in todo]) if s in res]
    report(a, out, order, res, prov, genomes, baitname, bait, baitinfo, start, status, todo)


def ev(x):
    """BLAST reports 0.0 when the E-value underflows; printing '0' would overstate it."""
    return "<1e-180" if x <= 0 else "%.1g" % x


def report(a, out, order, res, prov, genomes, baitname, bait, baitinfo, start, status, todo):
    def row(c):
        return ["%.0f" % c["score"], "%.1f" % c["pid"], "%.0f%%" % (100 * c["cov"]), ev(c["evalue"]),
                len(c["hsps"]), c["stops"], c["gene"], "%s (%s) %s:%d-%d" % (c["chrom"], c["strand"], c["acc"],
                                                                            c["lo"], c["hi"])] if c else ["-"] * 8

    cols = "bits\t%s identity %%\t%s coverage\t%s E-value\t%s pieces\t%s internal stops\t%s gene\t%s location"
    with open(os.path.join(out, "FIGURE5.tsv"), "w") as fh:
        fh.write("species\tassembly accession\tcall\t"
                 + "best AGTR2 model: " + cols % (("AGTR2",) * 7) + "\t"
                 + "best non-AGTR2 model: " + cols % (("other",) * 7)
                 + "\tmodels examined\tlowest model examined (bits)\n")
        for s in order:
            r = res[s]
            fh.write("\t".join(map(str, [s, T.accession(prov.get(s, "")), r["call"]]
                                   + row(r["target"]) + row(r["paralog"])
                                   + [r["nmodels"], "%.0f" % r["lowest"]])) + "\n")
        for g in todo:
            if g["sp"] not in res:
                fh.write("%s\t%s\tSKIPPED%s\n" % (g["sp"], T.accession(prov.get(g["sp"], "")),
                                                  "\t%s" % status.get(g["sp"], "") + "\t-" * 17))
    with open(os.path.join(out, "FIGURE5_all_models.tsv"), "w") as fh:
        fh.write("species\trank\tgene this model really is\tis AGTR2\tbits\tidentity %\tcoverage\tpieces\t"
                 "internal stops\tchromosome\tscaffold\tstart\tend\tstrand\n")
        for s in order:
            for c in sorted(res[s]["models"], key=lambda c: -c["score"]):
                fh.write("\t".join(map(str, [s, c["rank"], c["gene"], "yes" if c["is_target"] else "no",
                                             "%.0f" % c["score"], "%.1f" % c["pid"], "%.3f" % c["cov"],
                                             len(c["hsps"]), c["stops"], c["chrom"], c["acc"], c["lo"], c["hi"],
                                             c["strand"]])) + "\n")

    # ---------------------------------------------------------------- controls
    exp = {g["sp"]: g["expect"] for g in genomes}
    floor = max([res[s]["paralog"]["score"] for s in order if res[s]["paralog"]] or [0])
    L = ["CONTROLS"]
    L.append("  PRESENT (AGTR2 reported in the literature; the screen must call it INTACT)")
    n_ok = n = 0
    for s in order:
        if exp.get(s) == "present":
            n += 1
            ok = res[s]["call"] == "INTACT"
            n_ok += ok
            L.append("      %-28s %-11s %s" % (s, res[s]["call"], "PASS" if ok else "FAIL"))
    if not n:
        L.append("      none in this run")
    L.append("  ABSENT (published loss; the screen must call it NO_AGTR2)")
    m_ok = m = 0
    for s in order:
        if exp.get(s) == "absent":
            m += 1
            ok = res[s]["call"] == "NO_AGTR2"
            m_ok += ok
            L.append("      %-28s %-11s %s" % (s, res[s]["call"], "PASS" if ok else "FAIL"))
    if not m:
        L.append("      none in this run")
    L.append("  AGTR1 FOUND (AGTR1 exists in every vertebrate; if the screen cannot find it, the screen is")
    L.append("  not working in that genome and its AGTR2 result means nothing)")
    a_ok = 0
    for s in order:
        c = res[s]["agtr1"]
        a_ok += bool(c)
        L.append("      %-28s %s" % (s, ("%s %.0f bits, %.1f%% id at %s" % (PARALOG, c["score"], c["pid"],
                                                                           c["acc"]) if c else "NOT FOUND  <-- FAIL")))
    L.append("  MEASURED PARALOGUE FLOOR (highest bit score of any NON-AGTR2 model in this run)")
    L.append("      %.0f bits  - a bar at or below this is a paralogue, not AGTR2" % floor)
    L.append("  SEARCH DEPTH (a NO_AGTR2 call means 'not among the models examined', so how far down the")
    L.append("  screen looked has to be stated. Either the list ran out before the --max-models cap - in")
    L.append("  which case nothing was cut off - or the lowest model examined must score BELOW the")
    L.append("  paralogue floor, so the search did not stop while a hit could still have been AGTR2.)")
    d_ok = 0
    for s in order:
        r = res[s]
        full = r["nmodels"] >= a.max_models
        ok = (r["lowest"] < floor) or not full
        d_ok += ok
        L.append("      %-28s %2d model(s), lowest %3.0f bits  %s" % (
            s, r["nmodels"], r["lowest"],
            "not truncated (every model found was examined)" if not full else
            "below the %.0f-bit floor" % floor if ok else
            "ABOVE the %.0f-bit floor  <-- FAIL, raise --max-models" % floor))
    L.append("  CALLS")
    for s in order:
        r = res[s]
        t = r["target"]
        L.append("      %-28s %-11s %s" % (s, r["call"], (
            "AGTR2 %.0f bits, %.1f%% id, %.0f%% cov, %d piece(s), %d stop(s)" % (
                t["score"], t["pid"], 100 * t["cov"], len(t["hsps"]), t["stops"]) if t else
            "no AGTR2 anywhere; best hit is %s at %.0f bits (%.1f%% id)" % (
                r["paralog"]["gene"], r["paralog"]["score"], r["paralog"]["pid"]) if r["paralog"] else "no hit")))
    txt = "\n".join(L)
    print("\n" + txt)
    open(os.path.join(out, "FIGURE5_controls.txt"), "w").write(txt + "\n")

    # ---------------------------------------------------------------- methods + log
    tver = T.run(["tblastn", "-version"]).splitlines()[0]
    npseudo = [s for s in order if res[s]["call"] == "PSEUDOGENE"]
    meth = (
        "METHODS (whole-genome AGTR2 orthology screen)\nGenerated by FIGURE5.py from this run.\n\n"
        "Bait. %s: %s\n\n"
        "Search. The bait was searched against each complete genome assembly with tblastn (%s; E < %s). "
        "High-scoring segment pairs were chained into collinear gene models on a single scaffold and strand, "
        "allowing introns of up to %s bp. Chaining is necessary because a pseudogene interrupted by stop codons "
        "and frameshifts is recovered in several pieces, none of which need outscore a single contiguous "
        "alignment to a paralogue; the summed score of the chain is what identifies it.\n\n"
        "Identification. The conceptual translation of each gene model was searched with blastp against the "
        "complete human RefSeq proteome, and the model was accepted as AGTR2 only if its best match was AGTR2 "
        "itself and no other human gene equalled or exceeded that score. Models matching another gene are "
        "reported under that gene; FIGURE5_all_models.tsv lists every "
        "model with its assignment. Each species was then called INTACT (an AGTR2 orthologue recovering at least "
        "%.0f%% of the bait in one piece with no internal stop codon), PSEUDOGENE (an AGTR2 orthologue carrying "
        "internal stop codons or recovered in several pieces)%s, or NO_AGTR2 (no AGTR2 orthologue anywhere in the "
        "assembly; the best hit is reported with the gene it belongs to).\n\n"
        "Paralogue floor. For every assembly the highest-scoring non-AGTR2 model is reported with its gene and "
        "bit score. The threshold separating an AGTR2 orthologue from a paralogue that merely resembles one is "
        "the highest of these across the run (%.0f bits), so it is measured from the data rather than chosen.\n\n"
        "Search depth. An absence can only be reported to the depth the screen reached, so that depth is "
        "stated rather than assumed: the %d highest-scoring gene models were identified in each assembly, and "
        "the lowest model examined scored between %.0f and %.0f bits - in every assembly below the %.0f-bit "
        "paralogue floor, so the search did not stop while a hit could still have been called AGTR2.\n\n"
        "Controls. Species whose AGTR2 is reported in the literature were required to be called INTACT (%d of "
        "%d); species with a published loss were required to be called NO_AGTR2 (%d of %d). As a control on the "
        "search itself, AGTR1 - which is present in every vertebrate - was required to be found in every "
        "assembly (%d of %d).\n" % (
            baitname, baitinfo, tver, EVALUE, format(a.max_intron, ","), 100 * a.min_coverage,
            "" if not npseudo else " - %s in this run" % ", ".join(npseudo),
            floor, a.max_models, min(res[s]["lowest"] for s in order), max(res[s]["lowest"] for s in order),
            floor, n_ok, n, m_ok, m, a_ok, len(order)))
    open(os.path.join(out, "FIGURE5_methods.txt"), "w").write(meth)
    with open(os.path.join(out, "FIGURE5_run_log.txt"), "w") as fh:
        fh.write("FIGURE5.py md5 %s\nTABLE1.py md5 %s (%s)\nstarted %s finished %s\n%s\n" % (
            T.md5(os.path.abspath(__file__)), T.md5(_T1), _T1, start.isoformat(timespec="seconds"),
            datetime.datetime.now().isoformat(timespec="seconds"), tver))
        fh.write("command: %s\nbait: %s, %d aa\n\n" % (" ".join(sys.argv), baitname, len(bait)))
        for s in order:
            fh.write("%-28s %s\n" % (s, prov.get(s, "")))
        fh.write("\n" + txt + "\n")

    drawn = False
    if a.no_figure:
        print("\n(--no-figure: nothing drawn. The tables below are the data.)")
    else:
        try:
            figure(out, order, res, floor, baitname, a.divergent_below)
            drawn = True
        except ImportError:
            print("\n(matplotlib not installed - figure skipped; FIGURE5.tsv has the numbers)")
    print("\nResults in %s: FIGURE5.tsv, FIGURE5_all_models.tsv, FIGURE5_controls.txt, FIGURE5_methods.txt, "
          "FIGURE5_run_log.txt%s" % (out, ", FIGURE5_orthology_scan.png/.pdf" if drawn else ""))
    # fingerprint, so console output can never be mistaken for an earlier run's scrollback
    print("THIS RUN: FIGURE5.py md5 %s, finished %s" % (
        T.md5(os.path.abspath(__file__)), datetime.datetime.now().isoformat(timespec="seconds")))


DARKRED, ORANGE, GOLD, PALE = "#8B1A1A", "#FF9A1F", "#F2C230", "#F7D4BE"
COL = {"INTACT": ORANGE, "DIVERGENT": PALE, "PSEUDOGENE": GOLD, "NO_AGTR2": DARKRED}
LAB = {"INTACT": "Intact AGTR2", "DIVERGENT": "Intact AGTR2, < %.0f%% identity",
       "PSEUDOGENE": "AGTR2 pseudogene (stops / in pieces)",
       "NO_AGTR2": "No AGTR2 — bar shows the best paralogue"}
# Hatched so a bar that is NOT the gene cannot be read as faint AGTR2 signal at a glance.
HATCH = {"NO_AGTR2": "///"}


def shade(call, c, cut):
    """Which colour key a bar gets. Only ever a presentation split of INTACT; never changes a call."""
    return "DIVERGENT" if (call == "INTACT" and cut > 0 and c and c["pid"] < cut) else call


def figure(out, order, res, floor, baitname, cut=0.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    # No clade headings are drawn. A heading sits between two rows, so which row it belongs to is
    # ambiguous whenever the species order is changed with --species; the species names carry the
    # grouping already, and a label that can be read against the wrong row is worse than none.
    ys = list(range(len(order)))[::-1]
    vals = [(res[s]["target"] or res[s]["paralog"]) for s in order]
    fig, ax = plt.subplots(figsize=(12.5, 0.55 * len(order) + 2.4))
    xmax = max([c["score"] for c in vals if c] + [floor]) * 1.75

    keys = [shade(res[s]["call"], c, cut) for s, c in zip(order, vals)]
    for y, s, c, k in zip(ys, order, vals, keys):
        ax.barh(y, c["score"] if c else 0, height=0.62, color=COL[k], hatch=HATCH.get(k),
                edgecolor="black", linewidth=0.6, zorder=3)
        if c:
            note = "Identity: %.1f%% | Coverage: %.1f%% | E-val: %s" % (
                c["pid"], 100 * c["cov"], ev(c["evalue"]) if c["evalue"] <= 0 else "%.2e" % c["evalue"])
            if res[s]["call"] == "PSEUDOGENE":
                note += " | %d piece(s), %d stop(s)" % (len(c["hsps"]), c["stops"])
            if res[s]["call"] == "NO_AGTR2":
                note += " | this hit is %s" % c["gene"]
            ax.text(c["score"] + xmax * 0.012, y, note, va="center", fontsize=8.5,
                    fontweight="bold", color=DARKRED if k == "DIVERGENT" else COL[k])
    ax.axvline(floor, ls="--", color="#555555", lw=1.3, zorder=2)

    ax.set_yticks(ys)
    ax.set_yticklabels([s.replace("_", " ") for s in order], fontsize=10, fontstyle="italic")
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel("BLAST+ bit score of the best chained model", fontsize=11, fontweight="bold")
    ax.set_title("AGTR2 Genomic Orthology Screening Profiles Across Target Taxa",
                 fontsize=13, fontweight="bold", color=DARKRED, pad=14)
    ax.grid(axis="x", ls="--", alpha=0.35, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    used = [k for k in ("INTACT", "DIVERGENT", "PSEUDOGENE", "NO_AGTR2") if k in keys]
    ax.legend([Patch(facecolor=COL[k], edgecolor="black", hatch=HATCH.get(k)) for k in used]
              + [Line2D([], [], ls="--", color="#555555")],
              [(LAB[k] % cut if k == "DIVERGENT" else LAB[k]) for k in used]
              + ["measured paralogue floor (%.0f bits)" % floor],
              loc="lower right", fontsize=9, frameon=True, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "FIGURE5_orthology_scan.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(out, "FIGURE5_orthology_scan.pdf"), bbox_inches="tight")


if __name__ == "__main__":
    main()
