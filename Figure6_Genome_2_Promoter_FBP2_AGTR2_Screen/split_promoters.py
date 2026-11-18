#!/usr/bin/env python3
"""
split_promoters.py - turn a plain text file of promoter sequences into one FASTA per promoter.

The sequences are read, not retyped. Each block in the input starts with a line naming the species and
the gene (for example "Myotis myotis FBP2 promoter:") and is followed by sequence lines. The script
normalises the name, writes one .fna per block, and prints the length and md5 of each so the files can be
checked against the source.

It refuses anything it cannot name confidently, rather than guessing - a promoter written to the wrong
filename is exactly the error this project has already had once.

    python split_promoters.py AGTR2_FBP2_Promoter.txt queries
"""
import hashlib
import os
import re
import sys

GENES = ("AGTR2", "FBP2")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    src, dst = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    os.makedirs(dst, exist_ok=True)

    blocks, cur = [], None
    for raw in open(src, encoding="utf-8", errors="replace"):
        ln = raw.strip()
        if not ln:
            continue
        # a header line names a gene and is not sequence; sequence lines are only ACGTN
        if re.fullmatch(r"[ACGTNacgtn]+", ln) and cur is not None:
            cur["seq"].append(ln.upper())
            continue
        gene = next((g for g in GENES if g.lower() in ln.lower()), None)
        if gene is None:
            if cur is not None and re.fullmatch(r"[ACGTNacgtn]+", ln):
                cur["seq"].append(ln.upper())
            continue
        words = re.sub(r"[^A-Za-z ]", " ", ln).split()
        sp = [w for w in words if w.lower() not in
              (gene.lower(), "promoter", "sequence", "region", "the", "of")]
        if len(sp) < 2:
            sys.exit("ERROR: cannot read a species name from: %r" % ln)
        species = "%s_%s" % (sp[0].capitalize(), sp[1].lower())
        cur = dict(gene=gene, species=species, seq=[], header=ln)
        blocks.append(cur)

    if not blocks:
        sys.exit("ERROR: no promoter blocks found in %s" % src)

    print("%-34s %-28s %7s  %s" % ("file", "from", "bp", "md5"))
    for b in blocks:
        seq = "".join(b["seq"])
        if len(seq) < 200:
            sys.exit("ERROR: %s %s came out %d bp - the file was not parsed as expected"
                     % (b["species"], b["gene"], len(seq)))
        name = "%s_promoter_%s" % (b["gene"], b["species"])
        p = os.path.join(dst, name + ".fna")
        with open(p, "w", newline="\n") as fh:
            fh.write(">%s %s %s | supplied by the user in %s\n" %
                     (name, b["species"], b["gene"], os.path.basename(src)))
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
        print("%-34s %-28s %7d  %s" % (name + ".fna", b["header"][:28], len(seq),
                                       hashlib.md5(seq.encode()).hexdigest()))
    print("\n%d promoter(s) written to %s" % (len(blocks), dst))


if __name__ == "__main__":
    main()
