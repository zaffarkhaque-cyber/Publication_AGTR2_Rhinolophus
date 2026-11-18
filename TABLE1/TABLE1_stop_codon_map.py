#!/usr/bin/env python3
"""
TABLE1_stop_codon_map.py - draws the stop codons and reading-frame changes in the AGTR2 copies of Aselliscus,
Hipposideros armiger and Hipposideros larvatus on the human AGTR2 protein (NP_000677.2, 363 aa), with
AGTR2-intact bats and human as reference tracks, and prints the same information as a table.

Input: the results folder of TABLE1.py (TABLE1_disruptions.tsv and
       TABLE1_alignment.fasta).
Usage: python TABLE1_stop_codon_map.py <results folder> [output name without extension]
Output: <name>.png (300 dpi), <name>.pdf, <name>.tsv (position-by-species table)
Needs: Python 3 with matplotlib.
"""
import csv
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

FOCUS = ["Aselliscus_stoliczkanus", "Hipposideros_armiger", "Hipposideros_larvatus"]
REFS = ["Homo_sapiens", "Myotis_myotis", "Megaderma_spasma"]
STOP_C, FS_C, SHARED_C = "#C00000", "#2F5597", "#FFE699"


def fasta(path):
    d, cur = {}, None
    for ln in open(path):
        ln = ln.strip()
        if ln.startswith(">"):
            cur = ln[1:].split()[0]
            d[cur] = ""
        elif cur:
            d[cur] += ln
    return d


def italic(sp):
    g, s = sp.split("_", 1)
    return "%s. %s" % (g[0], s.replace("_", " "))


def main(res, name):
    rows = list(csv.DictReader(open(os.path.join(res, "TABLE1_disruptions.tsv")), delimiter="\t"))
    poscol = [c for c in rows[0] if c.startswith("human AGTR2 position")][0]
    fa = fasta(os.path.join(res, "TABLE1_alignment.fasta"))
    human = next(v for k, v in fa.items() if k.startswith("Homo_sapiens_") and k.endswith("_query"))
    L = len(human)

    ev = {s: [] for s in FOCUS}
    for r in rows:
        if r["species"] in ev and r["event"] in ("STOP", "FRAME_CHANGE"):
            m = re.match(r"(\d+)(?:-(\d+))?$", r[poscol])
            if not m:
                continue
            a = int(m.group(1))
            b = int(m.group(2) or a)
            ev[r["species"]].append(dict(type=r["event"], a=a, b=b, codon=r["codon"], gpos=r["genome coordinate of codon"]))

    # positions shared by >= 2 of the three (stops exact; frame changes overlapping within 3 aa)
    shared = []
    allev = [(s, e) for s in FOCUS for e in ev[s]]
    for s, e in allev:
        sp = {s}
        for s2, e2 in allev:
            if s2 != s and e2["type"] == e["type"] and (
                    (e["type"] == "STOP" and e2["a"] == e["a"]) or
                    (e["type"] == "FRAME_CHANGE" and e2["a"] <= e["b"] + 3 and e["a"] <= e2["b"] + 3)):
                sp.add(s2)
        e["shared"] = sorted(sp)
        if len(sp) > 1:
            shared.append((e["a"], e["b"]))

    # ---------------------------------------------------------------- figure
    tracks = [s for s in REFS if s in fa or s == "Homo_sapiens"] + FOCUS
    STEP = 1.7                                   # vertical distance between tracks
    fig, ax = plt.subplots(figsize=(14, 1.15 * len(tracks) + 2.0))
    for a, b in shared:
        ax.axvspan(a - 1.5, b + 1.5, color=SHARED_C, zorder=0, lw=0)

    def runs(mask):
        out, start = [], None
        for p in range(1, L + 2):
            on = p <= L and mask[p - 1]
            if on and start is None:
                start = p
            elif not on and start is not None:
                out.append((start, p - 1))
                start = None
        return out

    def stagger(xs, gap=9):
        """label level (0, 1, 2) for each x so that labels closer than 'gap' aa do not overlap"""
        lv, last = [], {}
        for x in xs:
            k = 0
            while k in last and x - last[k] < gap:
                k += 1
            last[k] = x
            lv.append(k)
        return lv

    for i, s in enumerate(reversed(tracks)):
        y = i * STEP
        seq = human if s == "Homo_sapiens" else fa.get(s, "")
        have = [p <= len(seq) and seq[p - 1] != "-" for p in range(1, L + 1)]
        ax.plot([0.5, L + 0.5], [y, y], color="#D9D9D9", lw=2, zorder=1)
        for a, b in runs(have):
            ax.plot([a - 0.5, b + 0.5], [y, y], color="#7F7F7F", lw=8, solid_capstyle="butt", zorder=2)
        if s not in ev:
            ax.text(L + 6, y, "intact", va="center", fontsize=8, color="#548235")
            continue
        stops = sorted([e for e in ev[s] if e["type"] == "STOP"], key=lambda e: e["a"])
        fss = sorted([e for e in ev[s] if e["type"] != "STOP"], key=lambda e: e["a"])
        for e, k in zip(stops, stagger([e["a"] for e in stops])):
            ax.plot(e["a"], y + 0.3, marker="v", ms=9, color=STOP_C, zorder=4)
            ax.text(e["a"], y + 0.45 + 0.38 * k, "%d %s" % (e["a"], e["codon"]), ha="center", va="bottom",
                    fontsize=7, color=STOP_C, fontweight="bold" if len(e["shared"]) > 1 else "normal")
        for e, k in zip(fss, stagger([(e["a"] + e["b"]) / 2 for e in fss], gap=14)):
            mid = (e["a"] + e["b"]) / 2
            ax.plot([mid, mid], [y - 0.3, y + 0.3], color=FS_C, lw=3, zorder=4)
            ax.text(mid, y - 0.4 - 0.32 * k, "fs %s" % (e["a"] if e["a"] + 1 >= e["b"] else "%d-%d" % (e["a"], e["b"])),
                    ha="center", va="top", fontsize=7, color=FS_C,
                    fontweight="bold" if len(e["shared"]) > 1 else "normal")
    ax.set_yticks([i * STEP for i in range(len(tracks))])
    ax.set_yticklabels([italic(s) for s in reversed(tracks)], fontstyle="italic", fontsize=10)
    ax.set_ylim(-1.2, (len(tracks) - 1) * STEP + 1.4)
    ax.set_xlim(-2, L + 25)
    ax.set_xlabel("Position on human AGTR2 protein (NP_000677.2, amino acids 1-%d)" % L, fontsize=10)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(handles=[Line2D([], [], marker="v", ls="", color=STOP_C, ms=9, label="stop codon (position, codon)"),
                       Line2D([], [], color=FS_C, lw=3, label="reading-frame change (fs)"),
                       Line2D([], [], color=SHARED_C, lw=8, label="same position in >= 2 hipposiderids"),
                       Line2D([], [], color="#D9D9D9", lw=2, label="no aligned sequence")],
              loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(name + ".png", dpi=300)
    fig.savefig(name + ".pdf")

    # ---------------------------------------------------------------- table
    keys = sorted({(e["type"], e["a"], e["b"]) for s in FOCUS for e in ev[s]}, key=lambda k: (k[1], k[0]))
    hdr = ["event", "human position", "human residue"] + FOCUS + ["shared by"]
    out = []
    for typ, a, b in keys:
        row = [typ, str(a) if a == b else "%d-%d" % (a, b), human[a - 1]]
        n = 0
        for s in FOCUS:
            hit = [e for e in ev[s] if e["type"] == typ and e["a"] == a and e["b"] == b]
            seq = fa.get(s, "")
            if hit:
                row.append(hit[0]["codon"] if typ == "STOP" else "frameshift")
                n += 1
            else:
                row.append("no (%s)" % (seq[a - 1] if a <= len(seq) else "-"))
        row.append(n)
        out.append(row)
    with open(name + ".tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(hdr)
        w.writerows(out)
    widths = [13, 15, 14] + [max(26, len(s) + 2) for s in FOCUS] + [9]
    print("".join(h.ljust(w) for h, w in zip(hdr, widths)))
    for r in out:
        print("".join(str(v).ljust(w) for v, w in zip(r, widths)))
    print("\n'no (X)' = this species has residue X at that position (no stop there); 'no (-)' = no aligned sequence.")
    print("Saved %s.png, %s.pdf, %s.tsv" % (name, name, name))


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else os.path.join(sys.argv[1], "TABLE1_stop_codon_map"))
