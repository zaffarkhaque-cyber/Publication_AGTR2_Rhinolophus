# NOX5_Fragmentation_IL1F10_PPEF2

Loss and fragmentation of three genes in rhinolophoid bats, from whole-genome search.

**Independent loss of *IL1F10* in Rhinolophoidea, *PPEF2* in *Rhinolophus*, and fragmentation of
*NOX5* in Rhinolophidae and Hipposideridae.**

Three genes, three different phylogenetic depths. A single event cannot produce three distributions,
which is what makes the losses independent rather than merely co-occurring:

| gene | Rhinolophidae | Hipposideridae | Megadermatidae | other bats + human |
|---|---|---|---|---|
| *IL1F10* | absent 4/4 | absent 2/2 | absent 2/2 | intact 7/7 |
| *NOX5* | fragmentary 4/4 | fragmentary 2/2 | intact 2/2 | intact 7/7 |
| *PPEF2* | absent 4/4 | intact 2/2 | 1 fragmentary, 1 intact | intact 7/7 |

## What is here

```
code/            every script, including the self-tests
manifests/       the design file and the manifests parsed from it
queries/         one FASTA per query sequence
results/         result tables, the control report and the run log
results/raw/     unmodified tblastn output, one file per genome
figures/         figures as PNG (300 dpi) and PDF, with their legends
MATERIALS_AND_METHODS.md
MANIFEST.md5     a checksum for every file above
LICENSE
```

## Reproducing it

The search needs BLAST+ and the NCBI genome archives named in `manifests/genomes.tsv`, which are not
redistributed here. With those in place:

```
python code/parse_design.py <design file>     # manifests and queries from the design file
bash   code/run_interrogate.sh                # the search across all genomes
python code/locus_gap_check.py                # ambiguous bases at each locus
python code/make_controls.py                  # the control report
python code/make_figures.py                   # the figures
```

`python code/test_interrogate.py` runs the control tests and needs nothing but Python.

Everything downstream of the search can be regenerated from the tables alone, without repeating it:
the control report, the paper table and the figures all read `results/PROTEINS.tsv` and
`results/LOCUS_GAPS.tsv` and compute nothing that is not already in them.

## How to read the results

Read `results/CONTROLS.txt` before the tables. It reports whether each query was recovered from its
own genome, whether the control genes were recovered from every genome, and whether each clade
carrying a claim has a control of its own clade. If those fail, the tables are not interpretable.

`results/PAPER_TABLE.txt` states, for each gene, what the numbers support and what they do not.

## Scope

This is a genomic analysis. It measures gene loss. It contains no expression data, no challenge
experiment and no measured phenotype, so any physiological consequence of these losses is a
hypothesis the results make worth testing, not a finding reported here.

Generated 2026-11-19 by `code/make_package.py`. 87 files.
