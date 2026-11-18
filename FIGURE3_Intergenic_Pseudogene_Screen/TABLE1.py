#!/usr/bin/env python3
"""
TABLE1.py

Locates AGTR2 and its X-chromosome neighbours PLS3 and SLC6A14 (plus the paralog AGTR1) in genome assemblies,
and tests whether AGTR2 is absent from an intact PLS3-SLC6A14 interval.

Method, per genome:
  1. tblastn of the human RefSeq proteins (AGTR2, PLS3, SLC6A14, AGTR1) against the whole assembly (E < 1e-5).
  2. Hits are chained into collinear gene models (same scaffold and strand, introns <= --max-intron bp).
  3. Orthology: the protein of each gene model is searched (blastp) against the complete human RefSeq proteome.
     A model is an ORTHOLOG only if its best match is the query gene itself; otherwise it is reported under the
     gene it matches best (e.g. an AGTR2-query hit whose best match is AGTR1 is an AGTR1 locus, not AGTR2).
  4. Primary PLS3 and SLC6A14 copies: orthologous models with >= --min-coverage of the protein, preferring a pair
     on the same scaffold. Every other orthologous copy and every paralog hit is listed in TABLE1_gene_models.tsv.
  5. If PLS3 and SLC6A14 are on the same scaffold, the DNA between them is read from the assembly and its length,
     number of N's and longest N run are reported, together with whether an AGTR2 ortholog lies inside it.
  6. AGTR2 call: FOUND (orthologous model, >= --min-coverage, no internal stop codon), FOUND_DISRUPTED (same, but
     with >= 1 in-frame stop codon >= 10 aa from both ends of a hit), PARTIAL (orthologous, lower coverage) or
     NOT_FOUND (no orthologous model; the best AGTR2-query hit and the gene it belongs to are reported).
  7. AGTR2 disruptions: every stop codon (with its codon read from the genome) and every change of reading frame
     in each AGTR2 copy, on human AGTR2 coordinates, and which of them are shared between species (a disruption at
     the same position in several species most simply means one inactivation in their common ancestor).
     Negative control: species expected to have intact AGTR2 must show none.
  Controls (column 'expected_AGTR2' of the genome list): genomes marked 'present' must give FOUND and genomes
  marked 'absent' (published losses) must give NOT_FOUND. Each is reported PASS / FAIL.

Inputs (all downloadable from NCBI Datasets):
  --reference   human GCF_000001405.40 zip with 'Protein (FASTA)' and 'Annotation features (GFF)'
  --genomes     tab-separated file: species <TAB> group <TAB> expected_AGTR2 <TAB> path to the NCBI genome zip
                (zip must contain the genomic .fna). Lines starting with # are ignored.
Requires: Python >= 3.8 (standard library only) and NCBI BLAST+ (tblastn, blastp, makeblastdb, blastdbcmd).

Usage:
  python TABLE1.py --reference human.zip --genomes genomes.tsv --out results
         [--db-dir blast_dbs] [--threads 8] [--blast-bin DIR] [--max-intron 1000000] [--min-coverage 0.6]
Outputs (--out): TABLE1.tsv, TABLE1_gene_models.tsv, TABLE1_controls.txt, TABLE1_run_log.txt, TABLE1_queries.faa, raw_tblastn/.
"""
import argparse
import datetime
import glob
import hashlib
import io
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

GENES = ["AGTR2", "PLS3", "SLC6A14", "AGTR1"]
TARGET, FLANKS, PARALOG = "AGTR2", ("PLS3", "SLC6A14"), "AGTR1"
EVALUE = "1e-5"
MAX_INTERVAL = 20000000     # do not read intervals longer than this (bp)
STOP_EDGE_AA = 10           # a stop must be at least this many residues inside its exon hit to be reported
SHARE_WINDOW = 3            # two frame changes count as the same lesion if they fall within this many residues


# ============================================================================ utilities
def fix_path(p):
    """Windows path -> /mnt/<drive>/... when run under WSL/Linux."""
    if p and os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", p):
        p = "/mnt/" + p[0].lower() + p[2:].replace("\\", "/")
    return p


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("command failed: %s\n%s" % (" ".join(cmd), r.stderr.strip()))
    return r.stdout


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def accession(prov):
    """NCBI assembly accession from the provenance string (…/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna …)."""
    m = re.match(r"(GC[AF]_\d+\.\d+)", prov or "")
    return m.group(1) if m else "-"


def read_fasta(path):
    d, cur, buf = {}, None, []
    for ln in open(path):
        if ln.startswith(">"):
            if cur:
                d[cur] = "".join(buf)
            cur, buf = ln[1:].split()[0], []
        else:
            buf.append(ln.strip())
    if cur:
        d[cur] = "".join(buf)
    return d


# ============================================================================ human reference
def load_reference(zip_path, out):
    """protein.faa + GFF -> {protein: seq}, {protein: gene symbol}. Cached in out/."""
    cf, ci = os.path.join(out, "reference_proteins.faa"), os.path.join(out, "reference_symbols.tsv")
    tag = "%s\t%d" % (os.path.basename(zip_path), os.path.getsize(zip_path))
    if os.path.exists(cf) and os.path.exists(ci):
        lines = open(ci).read().splitlines()
        if lines and lines[0] == "#" + tag:
            return read_fasta(cf), dict(l.split("\t") for l in lines[1:])
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        pick = lambda xs: sorted(xs, key=lambda n: ("GCF_" not in n, n))
        faa = pick([n for n in names if n.endswith("protein.faa")])
        gff = pick([n for n in names if n.endswith(".gff") or n.endswith(".gff3")])
        if not faa or not gff:
            sys.exit("ERROR: the reference zip needs protein.faa and a .gff file")
        print("Reading %s ..." % faa[0])
        prot, cur, buf = {}, None, []
        for ln in io.TextIOWrapper(z.open(faa[0]), encoding="utf-8", errors="replace"):
            if ln.startswith(">"):
                if cur:
                    prot[cur] = "".join(buf)
                cur, buf = ln[1:].split()[0], []
            else:
                buf.append(ln.strip())
        if cur:
            prot[cur] = "".join(buf)
        print("Reading %s (protein -> gene; a few minutes) ..." % gff[0])
        sym = {}
        for ln in io.TextIOWrapper(z.open(gff[0]), encoding="utf-8", errors="replace"):
            if "\tCDS\t" not in ln:
                continue
            a = ln.rstrip("\n").split("\t")[8]
            p, g = re.search(r"protein_id=([^;]+)", a), re.search(r"(?:^|;)gene=([^;]+)", a)
            if p and g and p.group(1) not in sym:
                sym[p.group(1)] = g.group(1).upper()
    with open(cf, "w") as fh:
        for k in sorted(prot):
            fh.write(">%s\n%s\n" % (k, prot[k]))
    with open(ci, "w") as fh:
        fh.write("#" + tag + "\n" + "".join("%s\t%s\n" % kv for kv in sorted(sym.items())))
    return prot, sym


def best_protein(gene, prot, sym):
    """Curated NP_ before predicted XP_, then longest, then accession."""
    c = [p for p, s in sym.items() if s == gene and p in prot]
    return sorted(c, key=lambda p: (not p.startswith("NP_"), -len(prot[p]), p))[0] if c else None


def build_proteome_db(prot, dbdir):
    """Built in --db-dir, not --out: BLAST memory-maps its databases, which fails on a cloud-synced folder
    (OneDrive/Dropbox placeholders) with 'Database memory map file error'."""
    db = os.path.join(dbdir, "human_proteome_db")
    tag = "%s %d proteins" % (sys.platform, len(prot))   # Windows- and Linux-built BLAST databases are not interchangeable
    src = db + ".source.txt"
    if not (os.path.exists(src) and open(src).read().strip() == tag):
        for f in glob.glob(db + ".*"):
            os.remove(f)
    if not os.path.exists(db + ".pin") and not os.path.exists(db + ".00.pin"):
        tmp = db + ".faa"
        with open(tmp, "w") as fh:
            for k in sorted(prot):
                fh.write(">%s\n%s\n" % (k, prot[k]))
        run(["makeblastdb", "-in", tmp, "-dbtype", "prot", "-out", db, "-title", "human RefSeq proteome"])
        os.remove(tmp)
        open(src, "w").write(tag + "\n")
    return db


# ============================================================================ genomes
def build_db(species, zip_path, dbdir):
    """Nucleotide BLAST database from the genomic .fna in an NCBI zip; reused only if built from the same zip."""
    db = os.path.join(dbdir, species)
    src = db + ".source.txt"
    zsize = os.path.getsize(zip_path)
    if glob.glob(db + "*.nsq") or os.path.exists(db + ".nal"):
        prov = open(src).read().strip() if os.path.exists(src) else ""
        if "(from %s, %d bytes)" % (os.path.basename(zip_path), zsize) not in prov:
            raise RuntimeError("database %s exists but was not built from %s; delete it or use another --db-dir"
                               % (db, zip_path))
        return db, prov
    with zipfile.ZipFile(zip_path) as z:
        fna = [n for n in z.namelist() if n.endswith("_genomic.fna") and "cds_from" not in n and "rna_from" not in n]
        fna.sort(key=lambda n: (not os.path.basename(n).startswith("GCF_"), n))
        if not fna:
            raise RuntimeError("no *_genomic.fna in %s" % zip_path)
        tmp = db + "_tmp.fna"
        print("  unpacking %s ..." % os.path.basename(fna[0]))
        with z.open(fna[0]) as s, open(tmp, "wb") as d:
            shutil.copyfileobj(s, d, 16 << 20)
    print("  building BLAST database ...")
    try:
        run(["makeblastdb", "-in", tmp, "-dbtype", "nucl", "-parse_seqids", "-out", db, "-title", species])
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    prov = "%s (from %s, %d bytes)" % (os.path.basename(fna[0]), os.path.basename(zip_path), zsize)
    open(src, "w").write(prov + "\n")
    return db, prov


def tblastn(qfa, db, threads, cache):
    if os.path.exists(cache):
        return open(cache).read()
    out = run(["tblastn", "-query", qfa, "-db", db, "-evalue", EVALUE, "-outfmt",
               "6 qseqid sacc qstart qend sstart send evalue bitscore pident length qlen qseq sseq",
               "-max_target_seqs", "50", "-max_hsps", "500", "-num_threads", str(threads)])
    open(cache + ".tmp", "w").write(out)
    os.replace(cache + ".tmp", cache)
    return out


def parse_hsps(txt):
    by_q = {}
    for ln in txt.splitlines():
        f = ln.split("\t")
        if len(f) < 13:
            continue
        ss, se = int(f[4]), int(f[5])
        by_q.setdefault(f[0], []).append(dict(acc=f[1], qs=int(f[2]), qe=int(f[3]), s_lo=min(ss, se),
                                              s_hi=max(ss, se), strand="+" if se >= ss else "-",
                                              evalue=float(f[6]), bits=float(f[7]), pid=float(f[8]),
                                              alen=int(f[9]), qseq=f[11], sseq=f[12]))
    return by_q


def chain_hsps(hsps, max_intron, max_chains=10, q_overlap=15):
    """Best collinear chains (dynamic programming) per scaffold and strand, sorted by summed bit score."""
    groups = {}
    for h in hsps:
        groups.setdefault((h["acc"], h["strand"]), []).append(h)
    chains = []
    for (acc, strand), hs in sorted(groups.items()):
        hs = sorted(hs, key=lambda h: (h["s_lo"], h["qs"]))
        used = [False] * len(hs)
        for _ in range(max_chains):
            n = len(hs)
            best, prev = [0.0] * n, [-1] * n
            for i in range(n):
                if used[i]:
                    continue
                best[i] = hs[i]["bits"]
                for j in range(i):
                    if used[j] or best[j] <= 0:
                        continue
                    if hs[i]["s_lo"] - hs[j]["s_hi"] > max_intron or hs[i]["s_lo"] < hs[j]["s_hi"] - 10:
                        continue
                    ok = (hs[i]["qs"] >= hs[j]["qe"] - q_overlap) if strand == "+" else \
                         (hs[i]["qe"] <= hs[j]["qs"] + q_overlap)
                    if ok and best[j] + hs[i]["bits"] > best[i]:
                        best[i], prev[i] = best[j] + hs[i]["bits"], j
            cand = [i for i in range(n) if not used[i]]
            if not cand:
                break
            k = max(cand, key=lambda i: best[i])
            path = []
            while k >= 0:
                path.append(hs[k])
                used[k] = True
                k = prev[k]
            path.reverse()
            chains.append(dict(acc=acc, strand=strand, hsps=path, score=sum(h["bits"] for h in path),
                               lo=min(h["s_lo"] for h in path), hi=max(h["s_hi"] for h in path)))
    chains.sort(key=lambda c: (-c["score"], c["acc"], c["lo"]))
    return chains


def coverage(hsps, qlen):
    cov = set()
    for h in hsps:
        cov.update(range(h["qs"], h["qe"] + 1))
    return len(cov) / qlen


def blastp_best(records, db, threads, sym):
    """{id: protein} -> {id: (best gene, bits, second gene or None, bits)} against the human proteome."""
    if not records:
        return {}
    fd, path = tempfile.mkstemp(suffix=".faa")
    with os.fdopen(fd, "w") as fh:
        for k, v in records.items():
            fh.write(">%s\n%s\n" % (k, v))
    try:
        out = run(["blastp", "-query", path, "-db", db, "-evalue", "1e-3", "-outfmt", "6 qseqid sseqid bitscore",
                   "-max_target_seqs", "25", "-num_threads", str(threads)])
    finally:
        os.remove(path)
    best = {}
    for ln in out.splitlines():
        q, s, b = ln.split("\t")
        g = sym.get(s, "UNKNOWN:" + s)
        d = best.setdefault(q, {})
        d[g] = max(d.get(g, 0.0), float(b))
    res = {}
    for q, d in best.items():
        r = sorted(d.items(), key=lambda x: (-x[1], x[0]))
        res[q] = (r[0][0], r[0][1], r[1][0] if len(r) > 1 else None, r[1][1] if len(r) > 1 else 0.0)
    return res


def seq_info(db, acc, cache):
    """(chromosome label, sequence length) from the assembly's own FASTA title."""
    if acc not in cache:
        t = run(["blastdbcmd", "-db", db, "-entry", acc, "-outfmt", "%l|%t"]).strip().split("|", 1)
        m = re.search(r"\b(?:chromosome|linkage group)[:\s]+([A-Za-z0-9_]+)", t[1] if len(t) > 1 else "",
                      re.IGNORECASE)
        cache[acc] = ("chr " + m.group(1) if m else "unplaced scaffold", int(t[0]))
    return cache[acc]


def interval_ns(db, acc, lo, hi):
    if hi < lo:
        return 0, 0, 0
    s = "".join(run(["blastdbcmd", "-db", db, "-entry", acc, "-range", "%d-%d" % (lo, hi),
                     "-outfmt", "%s"]).split()).upper()
    runs = [len(m.group(0)) for m in re.finditer(r"N+", s)]
    return len(s), sum(runs), max(runs) if runs else 0


# ============================================================================ AGTR2 disruptions
def hit_columns(h):
    """[(human position or None if bat insertion, bat residue, genome coordinate of the codon's first base)]"""
    out, qp, k = [], h["qs"], 0
    for qc, sc in zip(h["qseq"], h["sseq"]):
        g = None
        if sc != "-":
            g = h["s_lo"] + 3 * k if h["strand"] == "+" else h["s_hi"] - 3 * k
            k += 1
        out.append((qp if qc != "-" else None, sc, g, k - 1))
        if qc != "-":
            qp += 1
    return out, k


def frame_of(h):
    return (h["strand"], (h["s_lo"] % 3) if h["strand"] == "+" else (h["s_hi"] % 3))


def codon_at(db, acc, g, strand):
    lo, hi = (g, g + 2) if strand == "+" else (g - 2, g)
    s = run(["blastdbcmd", "-db", db, "-entry", acc, "-range", "%d-%d" % (lo, hi),
             "-strand", "plus" if strand == "+" else "minus", "-outfmt", "%s"])
    return "".join(s.split()).upper()


def disruption_report(hsps, qlen, db):
    """Stop codons and reading-frame changes of one gene locus, on human protein coordinates.
    hsps: every orthologous tblastn hit at the locus (hits in different frames over the same DNA included).
    Each human position is assigned to the highest-scoring hit covering it; a change of reading frame between
    consecutive assigned positions is a FRAME_CHANGE (frameshift) between those two positions."""
    best = {}
    for h in sorted(hsps, key=lambda h: (-h["bits"], h["s_lo"])):
        cols, _ = hit_columns(h)
        for qp, sc, g, _k in cols:
            if qp is not None and sc != "-" and qp not in best:
                best[qp] = (frame_of(h), sc)
    proj = "".join(best[p][1] if p in best else "-" for p in range(1, qlen + 1))
    events, seen = [], set()
    for h in hsps:
        cols, naa = hit_columns(h)
        for qp, sc, g, k in cols:
            if sc == "*" and g not in seen:
                seen.add(g)
                edge = k < STOP_EDGE_AA or k > naa - STOP_EDGE_AA - 1
                pos = qp if qp is not None else "insertion near %d" % next(
                    (c[0] for c in cols[cols.index((qp, sc, g, k)):] if c[0] is not None), qlen)
                events.append(dict(type="STOP", pos=pos, gpos=g, codon=codon_at(db, h["acc"], g, h["strand"]),
                                   note="within 10 aa of a hit end (could be hit running past a junction)"
                                   if edge else "inside hit"))
    ps = sorted(best)
    for p1, p2 in zip(ps, ps[1:]):
        if best[p1][0] != best[p2][0]:
            events.append(dict(type="FRAME_CHANGE", pos="%d-%d" % (p1, p2), gpos="-", codon="-",
                               note="reading frame %s%d -> %s%d" % (best[p1][0][0], best[p1][0][1],
                                                                     best[p2][0][0], best[p2][0][1])))
    return dict(events=events, proj=proj)


def shared_events(disr, window=SHARE_WINDOW):
    """Group events of different species that fall at the same human position (stops exact; frame changes within
    +/- window aa). Returns [(type, position label, [species])] for events seen in >= 1 species."""
    groups = []
    for sp, d in disr.items():
        for e in d["events"]:
            if e["type"] == "STOP" and isinstance(e["pos"], int):
                key = ("STOP", e["pos"], e["pos"])
            elif e["type"] == "FRAME_CHANGE":
                a, b = map(int, e["pos"].split("-"))
                key = ("FRAME_CHANGE", a, b)
            else:
                continue
            for g in groups:
                if g[0] == key[0] and (g[1] == key[1] if key[0] == "STOP" else
                                       (key[1] <= g[2] + window and g[1] <= key[2] + window)):
                    g[3].add(sp)
                    g[1], g[2] = min(g[1], key[1]), max(g[2], key[2])
                    break
            else:
                groups.append([key[0], key[1], key[2], {sp}])
    groups.sort(key=lambda g: (-len(g[3]), g[1]))
    return [(g[0], str(g[1]) if g[0] == "STOP" else "%d-%d" % (g[1], g[2]), sorted(g[3])) for g in groups]


# ============================================================================ one genome
def evaluate(sp, db, queries, qfa, qtag, pdb, sym, out, a):
    hs = parse_hsps(tblastn(qfa, db, a.threads, os.path.join(out, "raw_tblastn", "%s.%s.tsv" % (sp, qtag))))
    models, recs = [], {}
    for g in GENES:
        ch = chain_hsps(hs.get(g, []), a.max_intron)
        if ch:
            ch = [c for c in ch if c["score"] >= 0.2 * ch[0]["score"]][:10]
        for i, c in enumerate(ch):
            p = "".join(h["sseq"].replace("-", "").replace("*", "X") for h in sorted(c["hsps"], key=lambda h: h["qs"]))
            c.update(query=g, rank=i + 1, id="%s_%d" % (g, i + 1), cov=coverage(c["hsps"], len(queries[g])),
                     pid=sum(h["pid"] * h["alen"] for h in c["hsps"]) / sum(h["alen"] for h in c["hsps"]),
                     stops=sum(h["sseq"][10:-10].count("*") for h in c["hsps"]))
            if len(p) >= 25:
                recs[c["id"]] = p
            models.append(c)
    back = blastp_best(recs, pdb, a.threads, sym)
    info = {}
    for c in models:
        b = back.get(c["id"])
        c["best"] = b[0] if b else "no blastp match"
        c["second"] = b[2] if b else None
        c["ortholog"] = bool(b) and b[0] == c["query"] and not (b[2] and b[2] != c["query"] and b[3] >= b[1])
        c["chrom"], c["seqlen"] = seq_info(db, c["acc"], info)
    orth = {g: [c for c in models if c["query"] == g and c["ortholog"]] for g in GENES}
    good = {g: [c for c in orth[g] if c["cov"] >= a.min_coverage] for g in GENES}

    # primary flanks: same-scaffold pair with the highest summed score, else the best of each
    p3, s14 = good[FLANKS[0]], good[FLANKS[1]]
    pairs = [(x, y) for x in p3 for y in s14 if x["acc"] == y["acc"]]
    if pairs:
        fp = max(pairs, key=lambda xy: xy[0]["score"] + xy[1]["score"])
    else:
        fp = (p3[0] if p3 else None, s14[0] if s14 else None)
    r = dict(models=models, flank={FLANKS[0]: fp[0], FLANKS[1]: fp[1]})

    t = orth[TARGET]
    if t and t[0]["cov"] >= a.min_coverage:
        r["call"] = "FOUND_DISRUPTED" if t[0]["stops"] > 0 else "FOUND"
        r["agtr2"] = t[0]
    elif t:
        r["call"], r["agtr2"] = "PARTIAL", t[0]
    else:
        r["call"], r["agtr2"] = "NOT_FOUND", None
        anyhit = [c for c in models if c["query"] == TARGET]
        if anyhit:
            bh = anyhit[0]
            par = [c for c in orth[PARALOG] if c["acc"] == bh["acc"] and c["lo"] <= bh["hi"] and bh["lo"] <= c["hi"]]
            r["best_other"] = (bh, bool(par))
    if r["agtr2"]:
        p = r["agtr2"]
        locus = [h for c in t if c["acc"] == p["acc"] and c["strand"] == p["strand"]
                 and c["lo"] <= p["hi"] + 2000 and p["lo"] <= c["hi"] + 2000 for h in c["hsps"]]
        r["disr"] = disruption_report(locus, len(queries[TARGET]), db)

    # PLS3-SLC6A14 interval
    x, y = fp
    if x and y and x["acc"] == y["acc"]:
        a1, a2 = sorted([x, y], key=lambda c: c["lo"])
        lo, hi = a1["hi"] + 1, a2["lo"] - 1
        r["interval"] = dict(acc=x["acc"], lo=lo, hi=hi)
        if hi - lo + 1 > MAX_INTERVAL:
            r["interval"]["note"] = "longer than %d bp - not read" % MAX_INTERVAL
        else:
            r["interval"]["len"], r["interval"]["n"], r["interval"]["nrun"] = interval_ns(db, x["acc"], lo, hi)
        inside = [c for c in t if c["acc"] == x["acc"] and c["lo"] >= lo and c["hi"] <= hi]
        r["interval"]["agtr2_inside"] = "yes (%s:%d-%d)" % (inside[0]["acc"], inside[0]["lo"], inside[0]["hi"]) \
            if inside else "no"
    return r


def loc(c):
    if not c:
        return "-"
    return "%s (%s) %s:%d-%d" % (c["chrom"], c["strand"], c["acc"], c["lo"], c["hi"])


# ============================================================================ main
def main():
    ap = argparse.ArgumentParser(description="AGTR2 locus evaluation (Table 1)")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--db-dir", default=None)
    ap.add_argument("--proteome-db-dir", default=None,
                    help="where to build the human proteome BLAST database (default: --db-dir). Under WSL, point "
                         "this at a native Linux path such as ~/blastdb: BLAST memory-maps its databases and the "
                         "LMDB component is unreliable on a /mnt/c (9p) mount or a cloud-synced folder.")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--blast-bin", default=None, help="folder containing the BLAST+ programs, if not on the PATH")
    ap.add_argument("--max-intron", type=int, default=1000000)
    ap.add_argument("--min-coverage", type=float, default=0.60)
    a = ap.parse_args()
    start = datetime.datetime.now()
    if a.blast_bin:
        os.environ["PATH"] = fix_path(a.blast_bin) + os.pathsep + os.environ["PATH"]
    for t in ("tblastn", "blastp", "makeblastdb", "blastdbcmd"):
        if shutil.which(t) is None:
            sys.exit("ERROR: %s not found (install NCBI BLAST+ or give --blast-bin)" % t)
    out = fix_path(a.out)
    dbdir = fix_path(a.db_dir) if a.db_dir else os.path.join(out, "blast_dbs")
    for d in (out, dbdir, os.path.join(out, "raw_tblastn")):
        os.makedirs(d, exist_ok=True)
    ref = fix_path(a.reference)

    prot, sym = load_reference(ref, out)
    queries, qacc = {}, {}
    for g in GENES:
        p = best_protein(g, prot, sym)
        if not p:
            sys.exit("ERROR: %s not found in the reference annotation" % g)
        queries[g], qacc[g] = prot[p], p
    qfa = os.path.join(out, "TABLE1_queries.faa")
    with open(qfa, "w") as fh:
        for g in GENES:
            fh.write(">%s\n%s\n" % (g, queries[g]))
    qtag = md5(qfa)[:10]
    pdbdir = os.path.expanduser(fix_path(a.proteome_db_dir)) if a.proteome_db_dir else dbdir
    os.makedirs(pdbdir, exist_ok=True)
    pdb = build_proteome_db(prot, pdbdir)
    print("Queries: " + ", ".join("%s %s (%d aa)" % (g, qacc[g], len(queries[g])) for g in GENES))

    genomes = []
    for ln in open(fix_path(a.genomes), encoding="utf-8"):
        if ln.strip() and not ln.startswith("#"):
            f = [x.strip() for x in ln.rstrip("\n").split("\t")]
            if len(f) < 4 or f[2] not in ("present", "absent", "test"):
                sys.exit("ERROR: bad line in %s (need species, group, present/absent/test, accession, zip):\n%s"
                         % (a.genomes, ln))
            # accession is optional and may sit before or after the path; a 4-column file still works
            acc, zp = "", f[3]
            if len(f) >= 5:
                acc, zp = (f[3], f[4]) if re.match(r"GC[AF]_\d+\.\d+$", f[3]) else (f[4], f[3])
            genomes.append(dict(sp=f[0], group=f[1], expect=f[2], acc=acc, zip=fix_path(zp)))

    res, prov = {}, {}
    for g in genomes:
        print("\n==== %s" % g["sp"])
        if not os.path.exists(g["zip"]):
            prov[g["sp"]] = "SKIPPED - file not found: %s" % g["zip"]
            print("  " + prov[g["sp"]])
            continue
        try:
            db, prov[g["sp"]] = build_db(g["sp"], g["zip"], dbdir)
        except (RuntimeError, zipfile.BadZipFile) as e:
            prov[g["sp"]] = "SKIPPED - %s" % str(e).splitlines()[0]
            print("  " + prov[g["sp"]])
            continue
        got = accession(prov[g["sp"]])
        if g["acc"] and got != g["acc"]:
            sys.exit("ERROR: %s is listed in %s as assembly %s, but the file %s contains %s.\nThe genome list and "
                     "the file must agree - fix one of them." % (g["sp"], os.path.basename(a.genomes), g["acc"],
                                                                 os.path.basename(g["zip"]), got))
        r = res[g["sp"]] = evaluate(g["sp"], db, queries, qfa, qtag, pdb, sym, out, a)
        iv = r.get("interval")
        print("  AGTR2 %-9s %s" % (r["call"], loc(r["agtr2"])))
        for fl in FLANKS:
            print("  %-7s %s" % (fl, loc(r["flank"][fl])))
        if iv:
            print("  interval %s" % (iv.get("note") or "%d bp, %d N (longest run %d), AGTR2 inside: %s" % (
                iv["len"], iv["n"], iv["nrun"], iv["agtr2_inside"])))
        if "best_other" in r:
            bh, par = r["best_other"]
            print("  best AGTR2-query hit belongs to %s (%.1f%% id) at %s%s" % (
                bh["best"], bh["pid"], loc(bh), " = the AGTR1 ortholog locus" if par else ""))

    # ---------------------------------------------------------------- outputs
    with open(os.path.join(out, "TABLE1.tsv"), "w") as fh:
        fh.write("species\tassembly accession\tgroup\tAGTR2 call\tAGTR2 location\tAGTR2 coverage\tAGTR2 identity %\t"
                 "AGTR2 internal stop codons\tAGTR2 separate tblastn hits (1 = continuous reading frame)\t"
                 "PLS3 location\tSLC6A14 location\tPLS3 and SLC6A14 on same scaffold\tinterval bp\tN in interval\t"
                 "longest N run\tAGTR2 inside interval\tif AGTR2 not found: best AGTR2-query hit (gene, identity, "
                 "location)\tbest hit is the AGTR1 ortholog locus\tother orthologous PLS3 copies\t"
                 "other orthologous SLC6A14 copies\tassembly file\n")
        for g in genomes:
            s = g["sp"]
            if s not in res:
                fh.write("%s\t-\t%s\tSKIPPED%s\t%s\n" % (s, g["group"], "\t-" * 17, prov[s]))
                continue
            r = res[s]
            t, iv = r["agtr2"], r.get("interval", {})
            bo = r.get("best_other")
            extra = []
            for fl in FLANKS:
                prim = r["flank"][fl]
                extra.append("; ".join("%s cov %.0f%% %d hits" % (loc(c), 100 * c["cov"], len(c["hsps"]))
                                       for c in r["models"] if c["query"] == fl and c["ortholog"] and c is not prim)
                             or "-")
            fh.write("\t".join(map(str, [
                s, accession(prov[s]), g["group"], r["call"], loc(t), "%.0f%%" % (100 * t["cov"]) if t else "-",
                "%.1f" % t["pid"] if t else "-", t["stops"] if t else "-", len(t["hsps"]) if t else "-",
                loc(r["flank"][FLANKS[0]]), loc(r["flank"][FLANKS[1]]),
                "yes" if iv else "no", iv.get("len", iv.get("note", "-")), iv.get("n", "-"), iv.get("nrun", "-"),
                iv.get("agtr2_inside", "-"),
                "%s, %.1f%%, %s" % (bo[0]["best"], bo[0]["pid"], loc(bo[0])) if bo else "-",
                ("yes" if bo[1] else "no") if bo else "-", extra[0], extra[1], prov[s]])) + "\n")
    with open(os.path.join(out, "TABLE1_gene_models.tsv"), "w") as fh:
        fh.write("species\tassembly accession\tquery\tmodel rank\tchromosome\tscaffold\tstart\tend\tstrand\t"
                 "exon hits\tquery coverage\t"
                 "identity %\tsummed bits\tinternal stop codons\tbest match in human proteome\tsecond best\t"
                 "orthologous to query\n")
        for g in genomes:
            for c in res.get(g["sp"], {}).get("models", []):
                fh.write("\t".join(map(str, [g["sp"], accession(prov.get(g["sp"], "")), c["query"], c["rank"],
                                             c["chrom"], c["acc"], c["lo"], c["hi"],
                                             c["strand"], len(c["hsps"]), "%.3f" % c["cov"], "%.1f" % c["pid"],
                                             "%.0f" % c["score"], c["stops"], c["best"], c["second"] or "-",
                                             "yes" if c["ortholog"] else "no"])) + "\n")
    # ---------------------------------------------------------------- AGTR2 stop codons / frame changes
    disr = {g["sp"]: res[g["sp"]]["disr"] for g in genomes if g["sp"] in res and "disr" in res[g["sp"]]}
    with open(os.path.join(out, "TABLE1_disruptions.tsv"), "w") as fh:
        fh.write("species\tAGTR2 call\tevent\thuman AGTR2 position (%s)\tgenome coordinate of codon\tcodon\tnote\n"
                 % qacc[TARGET])
        for g in genomes:
            if g["sp"] in disr:
                ev = disr[g["sp"]]["events"]
                for e in ev:
                    fh.write("\t".join(map(str, [g["sp"], res[g["sp"]]["call"], e["type"], e["pos"], e["gpos"],
                                                 e["codon"], e["note"]])) + "\n")
                if not ev:
                    fh.write("%s\t%s\tnone\t-\t-\t-\tcontinuous reading frame, no stop codon\n"
                             % (g["sp"], res[g["sp"]]["call"]))
    with open(os.path.join(out, "TABLE1_alignment.fasta"), "w") as fh:
        fh.write(">Homo_sapiens_%s_query\n%s\n" % (qacc[TARGET], queries[TARGET]))
        for g in genomes:
            if g["sp"] in disr:
                fh.write(">%s %s\n%s\n" % (g["sp"], res[g["sp"]]["call"], disr[g["sp"]]["proj"]))
    sh = shared_events({s: d for s, d in disr.items() if d["events"]})
    S = ["AGTR2 DISRUPTIONS ON HUMAN AGTR2 COORDINATES (%s; stops: exact position; frame changes: within 3 aa)"
         % qacc[TARGET]]
    for typ, pos, sps in sh:
        S.append("  %-12s at %-9s in %d species: %s" % (typ, pos, len(sps), ", ".join(sps)))
    if not sh:
        S.append("  no stop codon or frame change in any AGTR2 copy")
    clean = [g["sp"] for g in genomes if g["expect"] == "present" and g["sp"] in disr]
    dirty = [s for s in clean if disr[s]["events"]]
    S.append("  NEGATIVE CONTROL (species expected to have intact AGTR2): %d of %d show no stop codon or frame "
             "change %s" % (len(clean) - len(dirty), len(clean), "PASS" if not dirty else "FAIL: " + ", ".join(dirty)))
    stxt = "\n".join(S)
    print("\n" + stxt)
    open(os.path.join(out, "TABLE1_shared_disruptions.txt"), "w").write(stxt + "\n")

    L = ["CONTROLS (expected_AGTR2 column of %s)" % os.path.basename(a.genomes)]
    for g in genomes:
        if g["expect"] == "test":
            continue
        r = res.get(g["sp"])
        if not r:
            L.append("  %-28s expected %-7s  NOT RUN (%s)" % (g["sp"], g["expect"], prov[g["sp"]]))
            continue
        ok = (r["call"] == "FOUND") if g["expect"] == "present" else (r["call"] == "NOT_FOUND")
        L.append("  %-28s expected %-7s  got %-9s  %s" % (g["sp"], g["expect"], r["call"], "PASS" if ok else "FAIL"))
    txt = "\n".join(L)
    print("\n" + txt)
    open(os.path.join(out, "TABLE1_controls.txt"), "w").write(txt + "\n")
    # ---------------------------------------------------------------- methods (generated from this run)
    tver = run(["tblastn", "-version"]).splitlines()[0]
    nan = [g["sp"] for g in genomes if g["sp"] in res]
    nsk = [g["sp"] for g in genomes if g["sp"] not in res]
    npass = sum(1 for g in genomes if g["expect"] != "test" and g["sp"] in res
                and ((res[g["sp"]]["call"] == "FOUND") if g["expect"] == "present"
                     else (res[g["sp"]]["call"] == "NOT_FOUND")))
    nctrl = sum(1 for g in genomes if g["expect"] != "test" and g["sp"] in res)
    meth = (
        "METHODS\n"
        "Generated by TABLE1.py from this run. Every value below is the value actually used; see "
        "TABLE1_run_log.txt for file checksums.\n\n"
        "Genome assemblies. %d genome assemblies were screened (NCBI accessions in TABLE1.tsv, column 'assembly "
        "accession')%s. No genome annotation was used.\n\n"
        "Gene identification. The human RefSeq proteins %s (assembly GCF_000001405.40) were searched against each "
        "assembly with tblastn (%s; E < %s). High-scoring segment pairs were chained into collinear gene models on "
        "a single scaffold and strand, allowing introns of up to %s bp. Orthology was assigned by reciprocal best "
        "hit: the conceptual translation of each gene model, and of each exon hit within it, was searched with "
        "blastp against the complete human RefSeq proteome, and a model was accepted as the orthologue only if its "
        "best match was the query gene itself; otherwise it is reported under the gene it matched (TABLE1_gene_"
        "models.tsv lists every model and its assignment). Models recovering at least %.0f%% of the query protein "
        "were called present.\n\n"
        "Absence. AGTR2 was called absent when no orthologous model was recovered. For each assembly the PLS3-"
        "SLC6A14 interval was located, and its length, total N content and longest run of N were measured directly "
        "from the assembly (blastdbcmd). Where AGTR2 was absent, the best remaining AGTR2-query hit is reported "
        "with the gene its translation matched and whether it coincides with the AGTR1 orthologue locus.\n\n"
        "Disruption. In-frame stop codons lying at least %d residues from both ends of an exon hit, and changes of "
        "reading frame between consecutively aligned residues, were recorded on human AGTR2 coordinates, with each "
        "stop codon read back from the genome. Lesions were counted as shared between species when stop codons fell "
        "at the identical residue and frame changes within %d residues.\n\n"
        "Controls. Species with an AGTR2 orthologue reported in the literature were required to return an intact "
        "gene, and species with published loss to return none: %d of %d passed. As a negative control for "
        "disruption calling, every species with an intact AGTR2 was required to show no stop codon and no frame "
        "change: %d of %d passed.\n" % (
            len(nan), "" if not nsk else "; %d further assemblies were unavailable and were skipped (%s)" % (
                len(nsk), ", ".join(nsk)),
            ", ".join("%s (%s)" % (g, qacc[g]) for g in GENES), tver, EVALUE, format(a.max_intron, ","),
            100 * a.min_coverage, STOP_EDGE_AA, SHARE_WINDOW, npass, nctrl, len(clean) - len(dirty), len(clean)))
    open(os.path.join(out, "TABLE1_methods.txt"), "w").write(meth)

    with open(os.path.join(out, "TABLE1_run_log.txt"), "w") as fh:
        fh.write("TABLE1.py md5 %s\nstarted %s finished %s\nPython %s, %s\n%s\n"
                 % (md5(os.path.abspath(__file__)), start.isoformat(timespec="seconds"),
                    datetime.datetime.now().isoformat(timespec="seconds"), sys.version.split()[0],
                    platform.platform(), run(["tblastn", "-version"]).splitlines()[0]))
        fh.write("command: %s\nreference %s md5 %s\ngenome list %s md5 %s\n" % (
            " ".join(sys.argv), ref, md5(ref), a.genomes, md5(fix_path(a.genomes))))
        fh.write("parameters: tblastn E < %s; max intron %d bp; min coverage %.2f\nqueries: %s\n\n" % (
            EVALUE, a.max_intron, a.min_coverage, ", ".join("%s %s" % (g, qacc[g]) for g in GENES)))
        for g in genomes:
            fh.write("%-28s %s\n" % (g["sp"], prov[g["sp"]]))
        fh.write("\n" + txt + "\n\n" + stxt + "\n")
    print("\nResults in %s: TABLE1.tsv, TABLE1_methods.txt, TABLE1_gene_models.tsv, TABLE1_controls.txt, "
          "TABLE1_disruptions.tsv, "
          "TABLE1_shared_disruptions.txt, TABLE1_alignment.fasta, TABLE1_run_log.txt" % out)


if __name__ == "__main__":
    main()
