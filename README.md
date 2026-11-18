# Publication_AGTR2_Rhinolophus

Code, inputs and published output for the computational analyses of **AGTR2 loss in *Rhinolophus***.

Every analysis here searches genome assemblies directly with BLAST+. **No genome annotation is used**,
so unannotated and newly released assemblies can be included, and no result depends on somebody else's
gene model being correct. Each analysis states its own controls and prints them *before* its results.

## What is in this repository

| Directory | Question it answers | Main output |
|---|---|---|
| [`TABLE1/`](TABLE1/) | Where are AGTR2 and its X-chromosome neighbours PLS3 and SLC6A14 in each assembly, and is AGTR2 missing from an otherwise intact PLS3–SLC6A14 interval? | `TABLE1.tsv`, stop-codon map |
| [`FIGURE3_Intergenic_Pseudogene_Screen/`](FIGURE3_Intergenic_Pseudogene_Screen/) | Where AGTR2 is not found, is a decayed remnant left between PLS3 and SLC6A14, or is the locus deleted outright? | `INTERGENIC.tsv`, `INTERGENIC_scan.png` |
| [`Figure5_AGTR2_Genomic_Orthology/`](Figure5_AGTR2_Genomic_Orthology/) | Is AGTR2 present **anywhere** in the genome — not just at its usual locus? | `FIGURE5.tsv`, `FIGURE5_orthology_scan.png` |

Each directory is self-contained: its own `README.md`, its own code, its own `results/` holding the
output of the published run, and its own generated Methods text. Start with the README inside the
directory you care about.

## How to read a result

Three conventions run through all three analyses, and they are the reason the calls can be checked
rather than taken on trust.

**Orthology is assigned by reciprocal best hit, not by score.** A candidate counts as AGTR2 only if,
when searched back against the *complete* human RefSeq proteome, its best match is AGTR2 and nothing
else ties or beats it. This is what separates AGTR2 from its paralogue AGTR1 — the two receptors are
only about one third identical, so a strong, highly significant BLAST hit can still be the wrong gene.

**Alignments are chained before they are judged.** A pseudogene broken by stop codons and frameshifts
comes back as several partial alignments, none of which need outscore a single contiguous alignment to
a paralogue. Methods that keep only the best single hit therefore report the paralogue and miss the
pseudogene. Collinear hits on one scaffold and strand are chained into one model first.

**Negative results state what they exclude.** An absence is only as good as the search behind it, so
each analysis reports how far down it looked — the number of models examined and the weakest score
reached, or a noise floor measured from shuffled sequence, and the assembly's N content across the
interval so that a gap is never read as a deletion.

## Reproducing a run

1. Download the human reference from NCBI Datasets, assembly **GCF_000001405.40**, selecting
   *Protein (FASTA)* and *Annotation features (GFF)*.
2. Download each assembly listed in that analysis's `genomes.tsv` (accession column) as
   *Genome sequences (FASTA)*. The convention in the packaged `genomes.tsv` is `genomes/<accession>.zip`;
   edit the last column to wherever you put them.
3. Edit the paths at the top of the `run_*.sh` script and run it.

Requirements: Python 3.8+ and NCBI BLAST+ 2.13+ on the `PATH` (`tblastn`, `blastp`, `makeblastdb`,
`blastdbcmd`); `matplotlib` for the figures.

Each `results/` directory holds the output of the published run. With the same assemblies and the same
command every file should be reproduced, apart from timestamps and paths in the run log. That log
records the **MD5 checksum of the script**, the BLAST+ version and every parameter, so what was run can
be checked against what is published here.

**BLAST databases must not sit in a cloud-synced folder** (OneDrive, Dropbox). BLAST memory-maps its
databases and placeholder files cannot be mapped; the symptom is `Database memory map file error`.

## Data availability

The genome assemblies are **not** redistributed here. Every one is public and is identified by its NCBI
accession in each analysis's `genomes.tsv`, which is also repeated in the results tables and the run
logs. Download them from NCBI under NCBI's terms.

## Licence

Code: MIT, see [LICENSE](LICENSE). Free to use, modify and redistribute, with attribution and without
warranty. The licence covers the code in this repository only, not the third-party genome assemblies
it analyses.

## Citation

> \<AUTHORS\>. \<TITLE\>. \<JOURNAL\> (\<YEAR\>). \<DOI\>

To cite the code itself rather than the paper, use the archived release:

> \<ZENODO DOI\>
