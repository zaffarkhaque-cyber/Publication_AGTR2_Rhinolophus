#!/usr/bin/env python3
"""
make_promoter_figure.py - draw the promoter figure from PROMOTER.tsv.

It reads the results table and nothing else, so the figure cannot disagree with the numbers: every bar
length, identity and call is taken straight from the file the analysis wrote. Re-running the screen is
not needed to redraw.

    python make_promoter_figure.py results
    python make_promoter_figure.py results --bat-agtr2 AGTR2_promoter_bat_recut

Two panels, as in the original experiment:
  left   FBP2  - the sensitivity control, searched across the whole genome
  right  AGTR2 - the question, searched at its own locus

Each species is measured against the query of its own class: bats against the bat-derived promoter,
birds against the bird-derived one. A bar is the number of query bases aligned; the label is identity.
Species where nothing was detected are drawn hatched at zero so an absence cannot be read as a short bar.
"""
import argparse
import csv
import os
import sys

DARKRED, ORANGE, GOLD, GREY = "#8B1A1A", "#FF9A1F", "#F2C230", "#BBBBBB"
BIRDGROUPS = ("hummingbirds", "other birds")


def load(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description="Draw the promoter figure from PROMOTER.tsv")
    ap.add_argument("results", help="the results directory holding PROMOTER.tsv")
    ap.add_argument("--bat-agtr2", default="AGTR2_promoter_bat_myotis")
    ap.add_argument("--bird-agtr2", default="AGTR2_promoter_bird_supplied")
    ap.add_argument("--bat-fbp2", default="FBP2_promoter_bat")
    ap.add_argument("--bird-fbp2", default="FBP2_promoter_bird")
    ap.add_argument("--out", default=None)
    ap.add_argument("--title",
                    default="Promoter retention at the $\\it{AGTR2}$ and $\\it{FBP2}$ loci "
                            "across bats and birds",
                    help="figure title; pass an empty string for none")
    a = ap.parse_args()

    rdir = os.path.abspath(a.results)
    rows = load(os.path.join(rdir, "PROMOTER.tsv"))
    if not rows:
        sys.exit("ERROR: PROMOTER.tsv is empty")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    # species order as the genome list gave it; birds last so the two classes read as blocks
    order, group = [], {}
    for r in rows:
        if r["species"] not in group:
            group[r["species"]] = (r["group"] or "").lower()
            order.append(r["species"])
    isbird = {s: any(b in group[s] for b in BIRDGROUPS) for s in order}
    order = [s for s in order if not isbird[s]] + [s for s in order if isbird[s]]

    by = {(r["species"], r["promoter"]): r for r in rows}

    def pick(sp, panel):
        q = (a.bird_fbp2 if isbird[sp] else a.bat_fbp2) if panel == "FBP2" else \
            (a.bird_agtr2 if isbird[sp] else a.bat_agtr2)
        return by.get((sp, q))

    missing = [sp for sp in order for p in ("FBP2", "AGTR2") if pick(sp, p) is None]
    if missing:
        sys.exit("ERROR: no row for %s - check the --*-agtr2 / --*-fbp2 names against PROMOTER.tsv"
                 % ", ".join(sorted(set(missing))))

    ys = list(range(len(order)))[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 0.42 * len(order) + 2.6), sharey=True)
    xmax = 1250

    for ax, panel, title in zip(axes, ("FBP2", "AGTR2"),
                                ("$\\it{FBP2}$ promoter — sensitivity control\n(searched across the whole genome)",
                                 "$\\it{AGTR2}$ promoter — at its own locus\n(upstream of $\\it{AGTR2}$; the PLS3–SLC6A14 interval where the gene is gone)")):
        for y, sp in zip(ys, order):
            r = pick(sp, panel)
            alen = int(r["aligned bp"] or 0)
            det = r["call"] == "DETECTED"
            col = (GOLD if isbird[sp] else ORANGE) if det else DARKRED
            ax.barh(y, alen if det else 28, height=0.62, color=col,
                    hatch=None if det else "///", edgecolor="black", linewidth=0.6, zorder=3)
            txt = ("%.1f%% ID" % float(r["identity %"])) if det else "not detected"
            ax.text((alen if det else 28) + 18, y, txt, va="center", fontsize=8.5,
                    fontweight="bold", color=col if det else DARKRED)
        ax.set_xlim(0, xmax)
        ax.set_ylim(-0.7, len(order) - 0.3)
        ax.set_xlabel("Query bases aligned", fontsize=10.5, fontweight="bold")
        ax.set_title(title, fontsize=11, fontweight="bold", color=DARKRED, pad=12)
        ax.grid(axis="x", ls="--", alpha=0.35, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    axes[0].set_yticks(ys)
    axes[0].set_yticklabels([s.replace("_", " ") for s in order], fontsize=9.5, fontstyle="italic")

    handles = [Patch(facecolor=ORANGE, edgecolor="black", label="detected, bat query"),
               Patch(facecolor=GOLD, edgecolor="black", label="detected, bird query"),
               Patch(facecolor=DARKRED, edgecolor="black", hatch="///",
                     label="not detected — nothing above the measured noise floor")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9.5,
               bbox_to_anchor=(0.5, -0.02))
    if a.title:
        fig.suptitle(a.title, fontsize=14, fontweight="bold", color=DARKRED, y=1.02)
    fig.tight_layout()

    stem = a.out or os.path.join(rdir, "PROMOTER_figure")
    fig.savefig(stem + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    print("Wrote %s.png and %s.pdf" % (stem, stem))
    print("FBP2 panel  : bats vs %s, birds vs %s" % (a.bat_fbp2, a.bird_fbp2))
    print("AGTR2 panel : bats vs %s, birds vs %s" % (a.bat_agtr2, a.bird_agtr2))


if __name__ == "__main__":
    main()
