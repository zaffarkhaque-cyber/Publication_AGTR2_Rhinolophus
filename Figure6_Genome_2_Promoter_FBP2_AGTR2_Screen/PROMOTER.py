#!/usr/bin/env python3
"""
PROMOTER.py - is the AGTR2 promoter still present, and how much of it, measured against chance?

This replaces a genome-wide blastn screen that took one HSP from one best hit. Four things were wrong
with that, and each is fixed here:

  1. NOTHING TIED THE HIT TO THE LOCUS. A blastn hit anywhere in a 2 Gb genome was read as "the promoter".
     Promoters are non-coding and repeat-rich, so that is not safe. Here the AGTR2 promoter is sought ONLY
     inside the PLS3-SLC6A14 interval - the place the gene belongs - and the FBP2 promoter only in the
     region immediately upstream of the FBP2 gene model. Both anchors are found in the genome itself by
     tblastn and reciprocal best hit, never by coordinates copied from another assembly.

  2. ONE HSP WAS TAKEN AS THE ANSWER. A decaying promoter survives in pieces. Collinear hits are chained
     and the TOTAL aligned length is reported, so decay is not mistaken for deletion.

  3. THERE WAS NO NOISE FLOOR, so a 28 bp hit at 100% identity was reported as a residual fragment. A 1 kb
     query against a 2 Gb genome produces such matches by chance constantly. Here each interval is shuffled
     20 times, preserving length and dinucleotide composition, and rescanned; the floor is the greater of
     the best shuffled score and mean + 3 SD. Anything at or below the floor is reported as NOT DETECTED,
     with the floor printed beside it.

  4. A NEGATIVE SAID NOTHING ABOUT WHAT IT EXCLUDED. The query is mutated from 5% to 60% divergence,
     embedded in shuffled sequence and put through the same scan, so every negative carries the divergence
     at which this method would still have found something.

The self-hit control is the one that matters most for this project: a query must recover ~100% of itself
over its full length in its OWN genome. A query labelled as one species but taken from another cannot do
that, and this control is what makes the discrepancy visible instead of silent.

    python PROMOTER.py --reference human.zip --genomes genomes.tsv --queries promoters.tsv --out results \\
                       --db-dir /path/to/blast_databases
"""
import argparse
import datetime
import json
import os
import random
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import TABLE1 as T                                                       # noqa: E402

SHUFFLES = 20            # shuffled replicates per interval, for the noise floor
GENOME_SHUFFLES = 10     # shuffled QUERIES for the floor of a whole-genome search (a genome cannot be
                         # shuffled; see genome_row). Fewer, because each one searches a whole genome.
FLOOR_SD = 3.0           # floor is max(shuffled) or mean + FLOOR_SD * SD, whichever is greater
DIVERGENCES = [5, 10, 20, 30, 40, 50, 60]      # % for the detection-limit curve
UPSTREAM = 10000         # bp searched upstream of an anchor gene for its promoter.
                         # tblastn finds the CODING start, not the transcription start, and the 5' UTR
                         # can be several kb. At 2000 the bird FBP2 promoter was not found even in Apus
                         # apus itself - a window too small to contain the thing being looked for.
EVALUE = "1e-3"
MAX_EVALUE = 1e-5        # a chained hit must also be significant in its own right
MIN_ALIGNED = 100        # below this a "remnant" is indistinguishable from a chance match in 300 kb.
                         # 50 was too loose: it admitted 77 bp hits in R. sinicus and R. affinis that
                         # only the mislabelled query found, and an 89 bp bat-to-bird hit in Apus.
SEED = 20261118


# --------------------------------------------------------------------------- sequence helpers
def get_seq(db, acc, lo, hi):
    """Plain uppercase sequence of one interval, straight from the assembly."""
    if hi < lo:
        return ""
    s = T.run(["blastdbcmd", "-db", db, "-entry", acc, "-range", "%d-%d" % (lo, hi), "-outfmt", "%s"])
    return "".join(s.split()).upper()


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def shuffle_preserving_dinuc(seq, rng):
    """Altschul-Erickson style shuffle: preserves length and dinucleotide composition, which ordinary
    base shuffling does not. Promoters are CpG- and repeat-biased, so a mononucleotide shuffle gives a
    floor that is far too low and makes chance hits look real."""
    s = [c for c in seq if c in "ACGT"]
    if len(s) < 10:
        return "".join(s)
    edges = {}
    for a, b in zip(s, s[1:]):
        edges.setdefault(a, []).append(b)
    for k in edges:
        rng.shuffle(edges[k])
    out, cur = [s[0]], s[0]
    for _ in range(len(s) - 1):
        if not edges.get(cur):
            cur = rng.choice([k for k in edges if edges[k]] or list(edges))
        out.append(edges[cur].pop())
        cur = out[-1]
    return "".join(out)


def mutate(seq, pct, rng):
    s = list(seq)
    for i in rng.sample(range(len(s)), int(len(s) * pct / 100.0)):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


# --------------------------------------------------------------------------- blastn against a sequence
def blastn_chain(query_path, subject_seq, task="blastn"):
    """Best collinear chain of blastn HSPs for one query against one interval.

    -subject is used instead of a database: the interval is a few hundred kb, the search is instant, and
    nothing is written to a shared database directory while other analyses read it.

    Two things here were wrong in the first version and are worth stating, because both produced
    confident nonsense rather than an error:

      * DUST (low-complexity filtering) was switched off. A promoter query containing a poly-A tract then
        matches poly-A runs scattered through every genome, and those hits were being reported as
        surviving promoter fragments. Filtering is ON.
      * Bit scores were summed over EVERY hit in the interval, so hundreds of hits to the same few query
        bases added up to thousands of bits for a 30 bp alignment. A promoter is one locus, so only ONE
        chain of collinear, co-located hits may contribute; the best such chain is reported.

    Returns (aligned query bp, weighted identity, chain bit score, best E-value, [HSPs of that chain])."""
    if not subject_seq:
        return 0, 0.0, 0.0, 1.0, []
    fd, sp = tempfile.mkstemp(suffix=".fna")
    with os.fdopen(fd, "w") as fh:
        fh.write(">interval\n%s\n" % subject_seq)
    try:
        out = subprocess.run(["blastn", "-query", query_path, "-subject", sp, "-task", task,
                              "-evalue", EVALUE, "-dust", "yes", "-soft_masking", "true",
                              "-outfmt", "6 qstart qend sstart send pident length bitscore evalue"],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout
    finally:
        os.remove(sp)
    hs = []
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 8:
            continue
        hs.append(dict(qs=int(f[0]), qe=int(f[1]), ss=int(f[2]), se=int(f[3]), pid=float(f[4]),
                       alen=int(f[5]), bits=float(f[6]), ev=float(f[7]),
                       strand="+" if int(f[3]) >= int(f[2]) else "-"))
    if not hs:
        return 0, 0.0, 0.0, 1.0, []

    return best_chain(hs)[:5]


def best_chain(hs):
    """The single best locus: hits grouped by sequence, strand and position, then a greedy collinear
    chain in which no query base is counted twice. A promoter may be broken into pieces; it may not be
    assembled out of fragments on different scaffolds or megabases apart.

    Returns (aligned query bp, weighted identity, bit score, best E-value, chain, sequence, start, end)."""
    span = 10 * max(h["qe"] for h in hs)
    best = None
    for acc, strand in sorted({(h.get("acc", ""), h["strand"]) for h in hs}):
        side = sorted([h for h in hs if h.get("acc", "") == acc and h["strand"] == strand],
                      key=lambda h: min(h["ss"], h["se"]))
        i = 0
        while i < len(side):
            grp, anchor = [side[i]], min(side[i]["ss"], side[i]["se"])
            j = i + 1
            while j < len(side) and min(side[j]["ss"], side[j]["se"]) - anchor <= span:
                grp.append(side[j])
                j += 1
            chain, usedq = [], set()
            for h in sorted(grp, key=lambda h: -h["bits"]):
                q = set(range(min(h["qs"], h["qe"]), max(h["qs"], h["qe"]) + 1))
                if q & usedq:
                    continue                     # the same query bases may not be counted twice
                chain.append(h)
                usedq |= q
            score = sum(h["bits"] for h in chain)
            if best is None or score > best[0]:
                best = (score, chain, usedq, acc,
                        min(min(h["ss"], h["se"]) for h in chain),
                        max(max(h["ss"], h["se"]) for h in chain))
            i = j
    score, chain, usedq, acc, lo, hi = best
    wid = sum(h["pid"] * h["alen"] for h in chain) / sum(h["alen"] for h in chain)
    return len(usedq), wid, score, min(h["ev"] for h in chain), chain, acc, lo, hi


def blastn_whole_genome(query_path, db, threads, task="blastn"):
    """The same search, against the ENTIRE genome rather than one interval.

    This is what a sensitivity control needs. FBP2 is not being asked "is it at its locus" - it is being
    asked "does this assay light up at all for a promoter that is still there". Anchoring that question
    to a window was wrong: it made the control fail in nine genomes including Apus apus against its own
    sequence, which is a statement about the window, not about the genome.

    Returns the same tuple as best_chain, or zeros when nothing is found."""
    out = subprocess.run(["blastn", "-query", query_path, "-db", db, "-task", task,
                          "-evalue", EVALUE, "-dust", "yes", "-soft_masking", "true",
                          "-num_threads", str(threads), "-max_target_seqs", "50",
                          "-outfmt", "6 sacc qstart qend sstart send pident length bitscore evalue"],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout
    hs = []
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 9:
            continue
        hs.append(dict(acc=f[0], qs=int(f[1]), qe=int(f[2]), ss=int(f[3]), se=int(f[4]),
                       pid=float(f[5]), alen=int(f[6]), bits=float(f[7]), ev=float(f[8]),
                       strand="+" if int(f[4]) >= int(f[3]) else "-"))
    if not hs:
        return 0, 0.0, 0.0, 1.0, [], "-", 0, 0
    return best_chain(hs)


def best_bits(query_path, subject_seq):
    return blastn_chain(query_path, subject_seq)[2]


def genome_row(sp, db, q, a, rng):
    """One whole-genome search, with a noise floor built by shuffling the QUERY.

    A 2 Gb genome cannot be shuffled, so the null is inverted: the query is shuffled, preserving length
    and dinucleotide composition, and each shuffle is searched against the same genome. That answers
    'what does a sequence of this length and composition score against this genome by chance', which is
    the right null for an unanchored search."""
    alen, pid, bits, ev, hsps, acc, lo, hi = blastn_whole_genome(q["path"], db, a.threads)
    sh = []
    for _ in range(a.genome_shuffles):
        fd, p = tempfile.mkstemp(suffix=".fna")
        with os.fdopen(fd, "w") as fh:
            fh.write(">shuffled\n%s\n" % shuffle_preserving_dinuc(q["seq"], rng))
        try:
            sh.append(blastn_whole_genome(p, db, a.threads)[2])
        finally:
            os.remove(p)
    mean = sum(sh) / len(sh) if sh else 0.0
    sd = (sum((x - mean) ** 2 for x in sh) / max(1, len(sh) - 1)) ** 0.5 if len(sh) > 1 else 0.0
    floor = max(max(sh) if sh else 0.0, mean + FLOOR_SD * sd)
    call = "DETECTED" if (bits > floor and ev < MAX_EVALUE and alen >= MIN_ALIGNED) else "NOT_DETECTED"
    return dict(sp=sp, q=q, call=call, note="whole genome, best locus %s:%d-%d" % (acc, lo, hi),
                span=0, ns=0, longest=0, alen=alen, pid=pid, bits=bits, ev=ev, floor=floor,
                acc=acc, lo=lo, hi=hi, hsps=hsps)


# --------------------------------------------------------------------------- anchors
def anchor_interval(sp, db, prot, sym, anchor, out, a, upstream=None):
    """Where to look, found in THIS genome rather than assumed from another assembly.

    anchor is either {"between": [geneA, geneB]} - the interval between two orthologues, used for AGTR2,
    whose own locus is deleted in the species of interest - or {"upstream_of": gene}, used for FBP2.
    Returns (accession, lo, hi, note) or (None, 0, 0, why not)."""
    # "upstream of the gene itself where the gene exists, the flanking interval where it does not".
    # A promoter belongs to its gene. Anchoring AGTR2's promoter to PLS3-SLC6A14 is only correct where
    # AGTR2 is gone; where the gene is present but has moved out of that interval - Hipposideros larvatus
    # is exactly this case - the flanking window searches the wrong place and returns a false absence.
    if "upstream_of_or_between" in anchor:
        spec = anchor["upstream_of_or_between"]
        acc, lo, hi, note = anchor_interval(sp, db, prot, sym, {"upstream_of": spec["gene"]},
                                            out, a, upstream)
        if acc is not None:
            return acc, lo, hi, note
        return anchor_interval(sp, db, prot, sym, {"between": spec["flanks"]}, out, a, upstream)

    genes = anchor.get("between") or [anchor["upstream_of"]]
    qfa = os.path.join(out, "anchors_%s.faa" % sp)
    with open(qfa, "w") as fh:
        for g in genes:
            p = T.best_protein(g, prot, sym)
            if not p:
                return None, 0, 0, "%s is not in the reference proteome" % g
            fh.write(">%s\n%s\n" % (g, prot[p]))
    # The cache key must name the GENES searched, not just the species. It did not, so the first anchor
    # searched for a species (PLS3+SLC6A14) was handed back for every later one, and FBP2 was reported
    # "not found" in all 14 genomes without ever being looked for.
    key = "%s.%s" % (sp, "-".join(sorted(genes)))
    hs = T.parse_hsps(T.tblastn(qfa, db, a.threads, os.path.join(out, "raw_tblastn", "%s.tsv" % key)))
    recs, chains = {}, {}
    for g in genes:
        ch = T.chain_hsps(hs.get(g, []), a.max_intron)
        if not ch:
            return None, 0, 0, "%s not found" % g
        chains[g] = ch[:5]
        for i, c in enumerate(ch[:5]):
            recs["%s_%d" % (g, i)] = "".join(h["sseq"].replace("-", "").replace("*", "X")
                                             for h in sorted(c["hsps"], key=lambda h: h["qs"]))
    back = T.blastp_best(recs, a.proteome_db, a.threads, sym)
    ok = {}
    for g in genes:
        for i, c in enumerate(chains[g]):
            b = back.get("%s_%d" % (g, i))
            if b and b[0] == g and not (b[2] and b[2] != g and b[3] >= b[1]):
                ok[g] = c
                break
        if g not in ok:
            return None, 0, 0, "%s found but reciprocal best hit names it %s" % (
                g, (back.get("%s_0" % g) or ["nothing"])[0])

    if "between" in anchor:
        x, y = ok[genes[0]], ok[genes[1]]
        if x["acc"] != y["acc"]:
            return None, 0, 0, "%s and %s are on different scaffolds" % tuple(genes)
        lo, hi = min(x["hi"], y["hi"]) + 1, max(x["lo"], y["lo"]) - 1
        if hi <= lo:
            return None, 0, 0, "flanking models overlap"
        return x["acc"], lo, hi, "between %s and %s" % tuple(genes)

    c = ok[genes[0]]
    up = upstream or a.upstream
    if c["strand"] == "+":
        lo, hi = max(1, c["lo"] - up), c["lo"] - 1
    else:
        lo, hi = c["hi"] + 1, c["hi"] + up
    return c["acc"], lo, hi, "%d bp upstream of %s (%s strand)" % (up, genes[0], c["strand"])


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="Promoter retention, measured against a noise floor")
    ap.add_argument("--reference", required=True, help="human NCBI zip: protein.faa + GFF")
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--queries", required=True,
                    help="TSV: label, search anchor json, source species, source anchor json, bp")
    ap.add_argument("--out", default="results")
    ap.add_argument("--db-dir", required=True)
    ap.add_argument("--proteome-db-dir", default=None)
    ap.add_argument("--species", default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--upstream", type=int, default=UPSTREAM)
    ap.add_argument("--max-intron", type=int, default=1000000)
    ap.add_argument("--shuffles", type=int, default=SHUFFLES)
    ap.add_argument("--genome-shuffles", type=int, default=GENOME_SHUFFLES,
                    help="shuffled queries used for the noise floor of a whole-genome search")
    ap.add_argument("--no-figure", action="store_true")
    a = ap.parse_args()

    out = T.fix_path(a.out)
    os.makedirs(os.path.join(out, "raw_tblastn"), exist_ok=True)
    dbdir = T.fix_path(a.db_dir)
    pdbdir = T.fix_path(a.proteome_db_dir) if a.proteome_db_dir else dbdir
    start = datetime.datetime.now()
    rng = random.Random(SEED)

    print("Reading the human reference ...")
    prot, sym = T.load_reference(T.fix_path(a.reference), out)
    a.proteome_db = T.build_proteome_db(prot, pdbdir)

    # queries, each with WHERE to look and WHERE IT CAME FROM. Provenance is required, not optional:
    # the species mix-up this screen is being rebuilt after was invisible precisely because the query
    # carried no record of the assembly and coordinates it was taken from.
    # Queries are not supplied as sequence. Each row says WHERE to take it from (a species and an anchor
    # gene) and WHERE to look for it (another anchor), and the program cuts it out of that genome itself.
    # Nothing is typed, so a query cannot be labelled as one species and taken from another - the mistake
    # this screen is being rebuilt after. The coordinates it used are recorded as the provenance.
    queries = []
    for ln in open(T.fix_path(a.queries), encoding="utf-8-sig"):
        if ln.lstrip("﻿").startswith("#") or not ln.strip():
            continue
        f = [x.strip() for x in ln.rstrip("\n").split("\t")]
        if len(f) < 5:
            sys.exit("ERROR: %s needs 5 columns: label, search anchor json, source species, "
                     "source anchor json, bp" % a.queries)
        # column 4 is either a FASTA you supply, or a json anchor telling the program where to cut it
        src4 = f[3].strip()
        queries.append(dict(label=f[0], anchor=json.loads(f[1]), src_sp=f[2],
                            src_anchor=json.loads(src4) if src4.startswith("{") else None,
                            src_file=None if src4.startswith("{") else T.fix_path(src4),
                            bp=int(f[4]) if len(f) > 4 and f[4].strip() else 0))

    # utf-8-sig, not utf-8: a file written by PowerShell's Set-Content -Encoding utf8 begins with a byte
    # order mark, which makes the first '#' comment look like data and fails three lines later with an
    # unhelpful IndexError. The line number is reported too, because "list index out of range" is not a
    # diagnosis.
    gl = []
    for i, ln in enumerate(open(T.fix_path(a.genomes), encoding="utf-8-sig"), 1):
        if ln.lstrip("﻿").startswith("#") or not ln.strip():
            continue
        f = [x.strip() for x in ln.rstrip("\n").split("\t")]
        if len(f) < 5:
            sys.exit("ERROR: %s line %d has %d column(s), expected 5 "
                     "(species, group, role, accession, path):\n  %r"
                     % (os.path.basename(a.genomes), i, len(f), ln.rstrip("\n")[:120]))
        gl.append(dict(sp=f[0], group=f[1], role=f[2], acc=f[3], zip=T.fix_path(f[4])))
    want = [s.strip() for s in a.species.split(",")] if a.species else None
    todo = [g for g in gl if want is None or g["sp"] in want]

    # ---------------------------------------------------------------- cut each query out of its genome
    os.makedirs(os.path.join(out, "queries"), exist_ok=True)
    for q in queries:
        if q["src_file"]:
            # You supplied the sequence. It is used exactly as given; nothing is trimmed or re-cut.
            # The self-hit control below is what checks it really came from the species it names.
            if not os.path.exists(q["src_file"]):
                sys.exit("ERROR: query file not found: %s" % q["src_file"])
            seqs = T.read_fasta(q["src_file"])
            if len(seqs) != 1:
                sys.exit("ERROR: %s must hold exactly one sequence (found %d)"
                         % (q["src_file"], len(seqs)))
            hdr, seq = list(seqs.items())[0]
            q.update(path=q["src_file"], seq=seq, header=hdr, md5=T.md5(q["src_file"]),
                     prov="supplied: %s (%d bp), claimed source %s"
                          % (os.path.basename(q["src_file"]), len(seq), q["src_sp"]))
            print("  query %-30s %5d bp  <- %s" % (q["label"], len(seq), q["prov"]))
            continue
        src = next((g for g in gl if g["sp"] == q["src_sp"]), None)
        if src is None:
            sys.exit("ERROR: source species %s for query %s is not in the genome list"
                     % (q["src_sp"], q["label"]))
        sdb, _ = T.build_db(src["sp"], src["zip"], dbdir)
        acc, lo, hi, note = anchor_interval(src["sp"], sdb, prot, sym, q["src_anchor"], out, a, q["bp"])
        if acc is None:
            sys.exit("ERROR: cannot take query %s from %s: %s" % (q["label"], q["src_sp"], note))
        seq = get_seq(sdb, acc, lo, hi)
        if len(seq) < 100:
            sys.exit("ERROR: query %s came out %d bp - too short to be a promoter region"
                     % (q["label"], len(seq)))
        q["path"] = os.path.join(out, "queries", "%s.fna" % q["label"])
        with open(q["path"], "w") as fh:
            fh.write(">%s %s %s:%d-%d (%s)\n%s\n" % (q["label"], src["sp"], acc, lo, hi, src["acc"], seq))
        q.update(seq=seq, header=q["label"], md5=T.md5(q["path"]),
                 prov="%s %s %s:%d-%d (%s)" % (src["sp"], src["acc"], acc, lo, hi, note))
        print("  query %-24s %5d bp  <- %s" % (q["label"], len(seq), q["prov"]))

    rows, limits, prov = [], [], {}
    for g in todo:
        sp = g["sp"]
        print("\n==== %s" % sp)
        db, src = T.build_db(sp, g["zip"], dbdir)
        prov[sp] = src
        for q in queries:
            if q["anchor"].get("genome"):
                r = genome_row(sp, db, q, a, rng)
                rows.append(r)
                print("  %-28s %-13s %4d bp, %.1f%% id, %6.0f bits (floor %3.0f), E %-8.1g [%s]"
                      % (q["label"], r["call"], r["alen"], r["pid"], r["bits"], r["floor"], r["ev"],
                         r["note"]))
                continue
            acc, lo, hi, note = anchor_interval(sp, db, prot, sym, q["anchor"], out, a)
            if acc is None:
                rows.append(dict(sp=sp, q=q, call="NO_ANCHOR", note=note, span=0, ns=0, longest=0,
                                 alen=0, pid=0.0, bits=0.0, floor=0.0, acc="-", lo=0, hi=0, hsps=[]))
                print("  %-26s NO ANCHOR: %s" % (q["label"], note))
                continue
            span, ns, longest = T.interval_ns(db, acc, lo, hi)
            seq = get_seq(db, acc, lo, hi)
            alen, pid, bits, ev, hsps = blastn_chain(q["path"], seq)

            # noise floor from the SAME interval, shuffled: composition and length held constant
            sh = [best_bits(q["path"], shuffle_preserving_dinuc(seq, rng)) for _ in range(a.shuffles)]
            mean = sum(sh) / len(sh)
            sd = (sum((x - mean) ** 2 for x in sh) / max(1, len(sh) - 1)) ** 0.5
            floor = max(max(sh), mean + FLOOR_SD * sd)

            # Three tests, all of which must pass. Above the floor alone is not enough: with a floor of 0
            # any chance match qualifies, which is how a 30 bp poly-A hit became a "surviving fragment".
            # A gap can hide an absence. It cannot unfind something that WAS found, so NOT_READABLE is
            # only ever applied to a negative: the first version voided Myotis myotis's own promoter at
            # 1831 bits because a 20 kb N-run sat elsewhere in a 589 kb interval.
            if bits > floor and ev < MAX_EVALUE and alen >= MIN_ALIGNED:
                call = "DETECTED"
            elif longest > 1000:
                call = "NOT_READABLE"
            else:
                call = "NOT_DETECTED"
            rows.append(dict(sp=sp, q=q, call=call, note=note, span=span, ns=ns, longest=longest,
                             alen=alen, pid=pid, bits=bits, ev=ev, floor=floor, acc=acc, lo=lo, hi=hi,
                             hsps=hsps))
            print("  %-28s %-13s %4d bp, %.1f%% id, %6.0f bits (floor %3.0f), E %-8.1g [%s N]"
                  % (q["label"], call, alen, pid, bits, floor, ev, longest))

        # detection limit: what would this genome's interval have revealed, had a diverged copy been there?
        for q in queries:
            r = [x for x in rows if x["sp"] == sp and x["q"] is q][0]
            if r["call"] == "NO_ANCHOR" or q["anchor"].get("genome"):
                continue      # a whole-genome search has no interval to plant a diverged copy into
            base = shuffle_preserving_dinuc(get_seq(db, r["acc"], r["lo"], r["hi"]), rng)
            worst = None
            for d in DIVERGENCES:
                planted = base[:len(base) // 2] + mutate(q["seq"], d, rng) + base[len(base) // 2:]
                if best_bits(q["path"], planted) > r["floor"]:
                    worst = d
            limits.append(dict(sp=sp, q=q["label"], limit=worst, floor=r["floor"]))

    report(out, todo, queries, rows, limits, prov, start, a)
    print("\nTHIS RUN: PROMOTER.py md5 %s, finished %s" % (
        T.md5(os.path.abspath(__file__)), datetime.datetime.now().isoformat(timespec="seconds")))


def report(out, todo, queries, rows, limits, prov, start, a):
    w = lambda n: open(os.path.join(out, n), "w", encoding="utf-8", newline="\n")
    role = {g["sp"]: g["role"] for g in todo}

    with w("PROMOTER.tsv") as fh:
        fh.write("\t".join(["species", "group", "role", "assembly", "promoter", "call",
                            "aligned bp", "query bp", "% of query", "identity %", "bit score",
                            "E-value", "noise floor (bits)", "pieces", "anchor", "scaffold",
                            "interval start",
                            "interval end", "interval bp", "N in interval", "longest N run"]) + "\n")
        for g in todo:
            for r in [x for x in rows if x["sp"] == g["sp"]]:
                q = r["q"]
                fh.write("\t".join(str(x) for x in [
                    g["sp"], g["group"], g["role"], g["acc"], q["label"], r["call"], r["alen"],
                    len(q["seq"]), "%.0f%%" % (100.0 * r["alen"] / len(q["seq"])), "%.1f" % r["pid"],
                    "%.0f" % r["bits"], "%.1g" % r.get("ev", 1.0), "%.0f" % r["floor"],
                    len(r["hsps"]), r["note"], r["acc"],
                    r["lo"], r["hi"], r["span"], r["ns"], r["longest"]]) + "\n")

    # every aligned block, in query and interval coordinates - what a promoter diagram is drawn from
    with w("PROMOTER_blocks.tsv") as fh:
        fh.write("\t".join(["species", "promoter", "block", "of", "query from", "query to",
                            "interval from", "interval to", "scaffold coordinate from",
                            "scaffold coordinate to", "identity %", "bits"]) + "\n")
        for r in rows:
            for i, h in enumerate(sorted(r["hsps"], key=lambda h: min(h["qs"], h["qe"])), 1):
                fh.write("\t".join(str(x) for x in [
                    r["sp"], r["q"]["label"], i, len(r["hsps"]), h["qs"], h["qe"], h["ss"], h["se"],
                    r["lo"] + min(h["ss"], h["se"]) - 1, r["lo"] + max(h["ss"], h["se"]) - 1,
                    "%.1f" % h["pid"], "%.0f" % h["bits"]]) + "\n")

    with w("PROMOTER_detection_limit.tsv") as fh:
        fh.write("species\tpromoter\tdetected up to % divergence\tnoise floor (bits)\n")
        for d in limits:
            fh.write("%s\t%s\t%s\t%.0f\n" % (d["sp"], d["q"],
                                             d["limit"] if d["limit"] is not None else "not even at 5%",
                                             d["floor"]))

    lines = ["CONTROLS", "=" * 100, "",
             "SELF-HIT (a query must recover itself, in full, in its OWN genome; anything less means the",
             "query did not come from the assembly its label claims - this is the control that catches a",
             "query taken from one species and labelled as another)"]
    selffail = 0
    for q in queries:
        hits = [r for r in rows if r["q"] is q and r["sp"] == q["src_sp"]]
        if not hits:
            lines.append("    %-26s its own genome was not analysed - self-hit NOT TESTED" % q["label"])
            continue
        r = hits[0]
        frac = 100.0 * r["alen"] / len(q["seq"])
        ok = frac >= 99.0 and r["pid"] >= 99.5
        selffail += not ok
        lines.append("    %-26s %s: %.0f%% of the query at %.1f%% identity   %s"
                     % (q["label"], r["sp"], frac, r["pid"], "PASS" if ok else
                        "FAIL - the query is not from this assembly"))
    # The other half of the control structure. The sensitivity control (FBP2, whole genome) shows the
    # assay CAN light up. This shows it does not light up for just anything: bats and birds are ~300 My
    # apart and non-coding sequence does not survive that, so every cross-class search must be negative.
    # Together they are what make a Rhinolophus zero mean something.
    grp = {g["sp"]: (g["group"] or "").lower() for g in todo}

    def klass(name):
        return "bird" if ("bird" in grp.get(name, "")) else "mammal"

    lines += ["", "CROSS-CLASS SPECIFICITY (a bat promoter must NOT be found in a bird genome, nor the",
              "reverse; any detection here means the assay returns a positive for sequence that cannot",
              "be present)"]
    cross = 0
    for q in queries:
        qk = klass(q["src_sp"])
        other = [g["sp"] for g in todo if klass(g["sp"]) != qk]
        hits = [r for r in rows if r["q"] is q and r["sp"] in other and r["call"] == "DETECTED"]
        cross += len(hits)
        lines.append("    %-30s %-7s query vs %2d opposite-class genome(s): %s"
                     % (q["label"], qk, len(other),
                        "all negative  PASS" if not hits else
                        "FAIL - detected in " + ", ".join(h["sp"] for h in hits)))
    lines.append("    %s" % ("every cross-class search negative" if not cross else
                             "%d cross-class detection(s): no negative in this run is trustworthy "
                             "until explained" % cross))

    lines += ["", "NOISE FLOOR (per interval, from %d shuffles preserving length and dinucleotide "
              "composition; a whole-genome row's floor comes from %d shuffled QUERIES instead)"
              % (a.shuffles, a.genome_shuffles)]
    for r in rows:
        if r["call"] != "NO_ANCHOR":
            lines.append("    %-26s %-22s floor %4.0f bits; observed %4.0f  %s"
                         % (r["sp"], r["q"]["label"], r["floor"], r["bits"], r["call"]))
    lines += ["", "DETECTION LIMIT (what a NOT_DETECTED excludes)"]
    for d in limits:
        lines.append("    %-26s %-22s detected to %s divergence"
                     % (d["sp"], d["q"], ("%d%%" % d["limit"]) if d["limit"] is not None else "<5%"))
    lines += ["", "ASSEMBLY GAPS (a promoter missing inside a gap is not a deletion)"]
    for r in rows:
        if r["call"] != "NO_ANCHOR":
            lines.append("    %-26s %-22s %d bp interval, %d N, longest run %d"
                         % (r["sp"], r["q"]["label"], r["span"], r["ns"], r["longest"]))
    lines += ["", "QUERY PROVENANCE", ""]
    for q in queries:
        lines.append("    %-26s %d bp, md5 %s, from %s (header in file: %s)"
                     % (q["label"], len(q["seq"]), q["md5"], q["prov"], q["header"]))
    lines += ["", "ASSEMBLY PROVENANCE", ""]
    for g in todo:
        lines.append("    %-26s %s  <- %s" % (g["sp"], g["acc"], prov.get(g["sp"], "?")))
    with w("PROMOTER_controls.txt") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))

    with w("PROMOTER_methods.txt") as fh:
        fh.write("\n".join([
            "METHODS (generated by PROMOTER.py from this run's parameters; nothing typed by hand)", "",
            "Each promoter was searched where the question required. The AGTR2 promoter was sought "
            "immediately upstream of the AGTR2 orthologue in species that retain the gene, and within the "
            "interval between the PLS3 and SLC6A14 orthologues in species that have lost it, since a "
            "deleted gene offers nothing to measure from. The FBP2 promoter was searched across the entire "
            "genome without any anchor: as a sensitivity control it is asked whether the assay detects a "
            "promoter at all, not whether one sits at a particular locus. Anchor genes were located in "
            "each assembly independently with tblastn, chained into collinear models, and accepted only "
            "when a blastp search of the model against the complete human RefSeq proteome returned that "
            "gene as the best match. Windows extended %d bp upstream of an anchor gene." % a.upstream, "",
            "Sequence was read directly from the assembly and searched with blastn (task blastn, E < %s, "
            "DUST low-complexity filtering and soft masking enabled). Hits were grouped by sequence, "
            "strand and position so that only one locus contributes, and chained greedily so that no query "
            "base is counted twice; the total number of query bases aligned and the length-weighted "
            "identity of that chain are reported." % EVALUE, "",
            "A noise floor was measured for every search. For an anchored search the interval was shuffled "
            "%d times preserving length and dinucleotide composition and rescanned. For a whole-genome "
            "search, where the subject cannot be shuffled, the query was instead shuffled %d times the "
            "same way and each shuffle searched against that genome. In both cases the floor is the "
            "greater of the highest shuffled score and the mean plus %.0f standard deviations. A promoter "
            "was called detected only if its chained score exceeded that floor AND its E-value was below "
            "%g AND at least %d query bases aligned; the floor alone is insufficient, since a floor of "
            "zero would admit any chance match. Where nothing was detected and the longest run of N "
            "exceeded 1000 bp the locus was reported as unreadable rather than as a deletion; a gap can "
            "conceal an absence but cannot invalidate a positive."
            % (a.shuffles, a.genome_shuffles, FLOOR_SD, MAX_EVALUE, MIN_ALIGNED), "",
            "Three controls are reported with every run. Each query was searched against the genome it is "
            "claimed to come from and must recover essentially all of itself at full identity, which "
            "identifies a sequence labelled as one species but taken from another. Every query was "
            "searched against all genomes of the opposite class, where roughly 300 million years of "
            "divergence leaves no alignable non-coding sequence, so any detection would indicate the "
            "assay returning a positive for sequence that cannot be present. The N content and longest N "
            "run of every interval searched are reported.", "",
            "The divergence at which this method would still detect a promoter was measured per interval by "
            "mutating the query at %s%% and embedding it in shuffled sequence of the same composition." %
            ", ".join(str(d) for d in DIVERGENCES), "",
            "Script PROMOTER.py md5 %s; TABLE1.py md5 %s." % (
                T.md5(os.path.abspath(__file__)), T.md5(os.path.join(HERE, "TABLE1.py"))),
        ]) + "\n")

    with w("PROMOTER_run_log.txt") as fh:
        fh.write("PROMOTER.py md5 %s\n" % T.md5(os.path.abspath(__file__)))
        fh.write("TABLE1.py md5 %s (%s)\n" % (T.md5(os.path.join(HERE, "TABLE1.py")),
                                              os.path.join(HERE, "TABLE1.py")))
        for q in queries:
            fh.write("query %s md5 %s  %s\n" % (q["label"], q["md5"], q["prov"]))
        fh.write("command: %s\n" % " ".join(sys.argv))
        fh.write("random seed %d, %d shuffles, floor = max(shuffled) or mean + %.0f SD\n"
                 % (SEED, a.shuffles, FLOOR_SD))
        fh.write("started %s, finished %s\n" % (start.isoformat(timespec="seconds"),
                                                datetime.datetime.now().isoformat(timespec="seconds")))
        for g in todo:
            fh.write("%-26s %s  %s\n" % (g["sp"], g["acc"], prov.get(g["sp"], "?")))
    print("\nResults in %s: PROMOTER.tsv, PROMOTER_blocks.tsv, PROMOTER_detection_limit.tsv, "
          "PROMOTER_controls.txt, PROMOTER_methods.txt, PROMOTER_run_log.txt" % out)


if __name__ == "__main__":
    main()
