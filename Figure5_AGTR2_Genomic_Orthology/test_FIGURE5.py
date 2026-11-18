"""
Logic test for FIGURE5.py.

Not a biological test. Three synthetic assemblies are planted so that the RIGHT ANSWER IS KNOWN, and the
program is required to produce it:

  mockA  AGTR1 + an intact AGTR2 coding sequence, contiguous, no stops          -> must be INTACT
  mockB  AGTR1 + the same AGTR2 coding sequence broken into 3 pieces on the      -> must be PSEUDOGENE
         minus strand, separated by 5 kb spacers, with 3 in-frame stop codons
  mockC  AGTR1 only                                                             -> must be NO_AGTR2

mockC is the test that matters: the planted AGTR1-like protein is 60% diverged from the AGTR2 bait, so the
bait DOES hit it. The program must identify that hit as AGTR1 by reciprocal blastp and refuse to call it
AGTR2. mockB is the other test that matters: the pseudogene is in 3 pieces, each of which scores LESS than
the single contiguous AGTR1 alignment in the same genome, so a method that took the best single HSP would
report AGTR1 and call AGTR2 absent. Only chaining recovers it.

The AGTR1 planted in all three is also what the program's own AGTR1 control checks, so the control must pass
in all three.

Real tblastn, makeblastdb and blastp are used. The genome layer and the human proteome are replaced.
"""
import importlib.util
import os
import random
import subprocess
import tempfile
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = HERE
WORK = os.path.join(tempfile.mkdtemp(prefix="FIGURE5_test_"), "work")
os.makedirs(WORK, exist_ok=True)
random.seed(4)

CODON = {"A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA", "E": "GAA", "G": "GGT",
         "H": "CAT", "I": "ATT", "L": "CTT", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCT", "S": "TCT",
         "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT", "X": "AAT", "*": "TAA"}
AA = "ACDEFGHIKLMNPQRSTVWY"


def rt(p):
    return "".join(CODON.get(c, "AAT") for c in p)


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def rnd(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def mutate(p, frac):
    l = list(p)
    for i in random.sample(range(len(l)), int(frac * len(l))):
        l[i] = random.choice(AA)
    return "".join(l)


def fasta(path, recs):
    with open(path, "w") as fh:
        for k, v in recs.items():
            fh.write(">%s\n%s\n" % (k, "\n".join(v[i:i + 70] for i in range(0, len(v), 70))))


# ------------------------------------------------------------------ sequences
bait = ""
name = None
for ln in open(os.path.join(SRC, "my_baits.faa")):
    if ln.startswith(">"):
        if bait:
            break
        name = ln[1:].strip()
    elif name:
        bait += ln.strip()
assert len(bait) > 300, "bait not read"
agtr1_prot = mutate(bait, 0.60)                    # the paralogue: diverged enough to be a different gene
print("bait %s, %d aa; planted paralogue %d aa" % (name, len(bait), len(agtr1_prot)))

agtr2_cds = rt(bait)
agtr1_cds = rt(agtr1_prot)

# mockB: a DECAYED pseudogene - 35% diverged from the bait, in four pieces on the minus strand,
# each separated by a 5 kb spacer, each carrying an in-frame stop codon. The divergence is deliberate:
# it makes every single piece score LESS than the one contiguous AGTR1 alignment in the same genome,
# which is the situation in which a best-single-HSP method reports AGTR1 and misses the pseudogene.
pseudo_cds = rt(mutate(bait, 0.35))
q = (len(pseudo_cds) // 12) * 3
pieces = [pseudo_cds[:q], pseudo_cds[q:2 * q], pseudo_cds[2 * q:3 * q], pseudo_cds[3 * q:]]
pieces = [p[:120] + "TAA" + p[123:] for p in pieces]        # one stop at codon 41, well inside each piece
assert all(len(p) > 240 for p in pieces)
broken = rc(pieces[0] + rnd(5000) + pieces[1] + rnd(5000) + pieces[2] + rnd(5000) + pieces[3])

GEN = {
    "mockA": rnd(60000) + agtr1_cds + rnd(60000) + agtr2_cds + rnd(60000),
    "mockB": rnd(60000) + agtr1_cds + rnd(60000) + broken + rnd(60000),
    "mockC": rnd(60000) + agtr1_cds + rnd(60000),
}
dbdir = os.path.join(WORK, "dbs")
os.makedirs(dbdir, exist_ok=True)
for sp, seq in GEN.items():
    fa = os.path.join(dbdir, sp + ".fna")
    fasta(fa, {"scaf_" + sp: seq})
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-parse_seqids",
                    "-out", os.path.join(dbdir, sp), "-title", sp], check=True,
                   stdout=subprocess.DEVNULL)

# a mini "human proteome": the two genes that matter, plus decoys
prot = {"NP_AGTR2": bait, "NP_AGTR1": agtr1_prot}
sym = {"NP_AGTR2": "AGTR2", "NP_AGTR1": "AGTR1"}
for g in ("BDKRB1", "CXCR3", "LPAR6", "APLNR"):
    prot["NP_" + g] = "".join(random.choice(AA) for _ in range(350))
    sym["NP_" + g] = g

gl = os.path.join(WORK, "genomes.tsv")
with open(gl, "w") as fh:
    fh.write("mockA\tOther bats\tpresent\tGCA_000000001.1\t%s\n" % os.path.join(dbdir, "mockA.fna"))
    fh.write("mockB\tHipposideridae\ttest\tGCA_000000002.1\t%s\n" % os.path.join(dbdir, "mockB.fna"))
    fh.write("mockC\tHummingbirds\tabsent\tGCA_000000003.1\t%s\n" % os.path.join(dbdir, "mockC.fna"))
ref = os.path.join(WORK, "dummy_reference.zip")
open(ref, "w").write("not read: load_reference is replaced")

# ------------------------------------------------------------------ run
spec = importlib.util.spec_from_file_location("FIGURE5", os.path.join(SRC, "FIGURE5.py"))
F = importlib.util.module_from_spec(spec)
sys.modules["FIGURE5"] = F
spec.loader.exec_module(F)
F.T.load_reference = lambda z, o: (prot, sym)
F.T.build_db = lambda sp, zp, dd: (os.path.join(dbdir, sp), "%s_genomic.fna (planted)" % sp)

out = os.path.join(WORK, "results")
sys.argv = ["FIGURE5.py", "--reference", ref, "--genomes", gl, "--out", out,
            "--db-dir", dbdir, "--baits", os.path.join(SRC, "my_baits.faa"),
            "--species", "mockA,mockB,mockC", "--threads", "4"]
F.main()

# ------------------------------------------------------------------ assertions
rows = {}
for i, ln in enumerate(open(os.path.join(out, "FIGURE5.tsv"))):
    f = ln.rstrip("\n").split("\t")
    if i:
        rows[f[0]] = f
ctl = open(os.path.join(out, "FIGURE5_controls.txt")).read()
meth = open(os.path.join(out, "FIGURE5_methods.txt")).read()
models = [l.split("\t") for l in open(os.path.join(out, "FIGURE5_all_models.tsv")).read().splitlines()[1:]]

fails = []


def check(ok, msg):
    print(("  PASS  " if ok else "  FAIL  ") + msg)
    if not ok:
        fails.append(msg)


print("\nASSERTIONS")
check(rows["mockA"][2] == "INTACT", "mockA (contiguous, no stops) is called INTACT, got %s" % rows["mockA"][2])
check(rows["mockB"][2] == "PSEUDOGENE", "mockB (3 pieces, 3 stops) is called PSEUDOGENE, got %s" % rows["mockB"][2])
check(rows["mockC"][2] == "NO_AGTR2", "mockC (paralogue only) is called NO_AGTR2, got %s" % rows["mockC"][2])

check(int(rows["mockB"][7]) >= 3, "mockB AGTR2 model is recovered in >= 3 of the 4 planted pieces, got %s"
      % rows["mockB"][7])
check(int(rows["mockB"][8]) >= int(rows["mockB"][7]) - 1,
      "mockB reports an internal stop in (nearly) every piece recovered: %s stop(s) in %s piece(s)"
      % (rows["mockB"][8], rows["mockB"][7]))
check(rows["mockA"][8] == "0", "mockA AGTR2 model reports 0 internal stops, got %s" % rows["mockA"][8])

check(rows["mockC"][3] == "-", "mockC reports no AGTR2 model at all, got %s bits" % rows["mockC"][3])
check(rows["mockC"][17] == "AGTR1", "mockC's best hit is named as AGTR1, got %s" % rows["mockC"][17])
check(rows["mockA"][9] == "AGTR2" and rows["mockB"][9] == "AGTR2",
      "mockA and mockB AGTR2 models are named AGTR2 (%s, %s)" % (rows["mockA"][9], rows["mockB"][9]))

# the point of chaining: in mockB, is the chained AGTR2 above the AGTR1 alignment, and would a single
# HSP have lost to it?
mb = [m for m in models if m[0] == "mockB"]
t = [m for m in mb if m[3] == "yes"]
p = [m for m in mb if m[3] == "no" and m[2] == "AGTR1"]
check(bool(t) and bool(p), "mockB has both an AGTR2 model and an AGTR1 model")
if t and p:
    chained, par = float(t[0][4]), float(p[0][4])
    per_piece = chained / int(t[0][7])
    check(chained > par, "chained AGTR2 (%.0f bits) beats AGTR1 (%.0f bits)" % (chained, par))
    check(per_piece < par, "each AGTR2 piece alone (~%.0f bits) would have LOST to AGTR1 (%.0f bits) - "
                           "chaining is what finds the pseudogene" % (per_piece, par))

check("1 of 1" in ctl.split("PRESENT")[0] or "PASS" in ctl.split("ABSENT")[0], "PRESENT control recorded")
check(ctl.count("NOT FOUND") == 0, "AGTR1 control: AGTR1 found in all three assemblies")
check("PSEUDOGENE" not in ctl.split("MEASURED")[0].split("ABSENT")[1].split("AGTR1 FOUND")[0],
      "no ABSENT-control species was called PSEUDOGENE")
check("FAIL" not in ctl, "no control FAILs")

floor = float(ctl.split("bits  - a bar at or below")[0].split()[-1])
check(abs(floor - max(float(m[4]) for m in models if m[3] == "no")) < 1,
      "the floor printed (%.0f) is the highest non-AGTR2 model in the run" % floor)

check("md5" in meth and os.path.basename("my_baits.faa") in meth, "Methods names the bait file and its md5")
check("chain" in meth.lower() and "paralogue" in meth.lower(), "Methods explains chaining and the floor")
check(str(len(bait)) in meth, "Methods states the bait length")
for s in ("mockA", "mockB", "mockC"):
    check(s in open(os.path.join(out, "FIGURE5_run_log.txt")).read(), "run log records %s" % s)
check(os.path.exists(os.path.join(out, "FIGURE5_orthology_scan.png")), "figure written")

# --divergent-below is presentation only: the same run with it on must produce identical tables
png = os.path.join(out, "FIGURE5_orthology_scan.png")
before = open(os.path.join(out, "FIGURE5.tsv"), "rb").read()
os.remove(png)
sys.argv += ["--divergent-below", "95"]          # mockA is 100% identical, mockB ~64%: splits them
F.main()
check(open(os.path.join(out, "FIGURE5.tsv"), "rb").read() == before,
      "--divergent-below changes no number in FIGURE5.tsv")
check(os.path.exists(png) and os.path.getsize(png) > 10000, "figure redrawn with the divergent shade")

print("\n%d assertion(s) failed" % len(fails))
for m in fails:
    print("  " + m)
sys.exit(1 if fails else 0)
