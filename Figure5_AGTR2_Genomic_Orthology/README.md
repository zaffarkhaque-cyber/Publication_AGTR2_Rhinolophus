# Figure 5 — whole-genome AGTR2 orthology screen

**Question:** searching the **entire assembly**, not just the AGTR2 locus, is there an AGTR2 orthologue
anywhere — intact, pseudogenised or translocated — or is the best hit in that genome merely a paralogue
such as AGTR1?

This is the complement to the other two analyses. `../TABLE1/` asks whether AGTR2 sits where it should,
between PLS3 and SLC6A14. `../FIGURE3_Intergenic_AGTR2_PSEUDOGENE_SCreen/` asks whether a decayed remnant
is left in that interval. This directory asks whether the gene is anywhere in the genome at all, so a
negative result cannot be explained by a translocation or by a mis-scaffolded locus.

**This directory is self-contained.** The analysis, the gene-finding library it uses, the genome list, the
bait, the run script and the output of the published run are all here. `TABLE1.py` is shared between the
analyses and the copy in `../TABLE1/` is the master — if you change it there, copy it here as well.

## What the analysis does

1. **Searches the whole assembly.** The AGTR2 bait protein is searched against each complete genome with
   `tblastn` (protein against translated DNA), which stays sensitive to a remnant long after nucleotide
   comparison fails.
2. **Chains the hits into gene models.** High-scoring segment pairs on one scaffold and strand are chained
   collinearly, allowing introns of up to 1,000,000 bp. **Chaining is the whole point of this script.** A
   pseudogene broken by stop codons and frameshifts is recovered in several pieces, and no single piece need
   outscore one contiguous alignment to a paralogue. A method that takes the single best HSP therefore
   reports the paralogue and calls the pseudogene absent. The summed score of the chain is what finds it.
3. **Identifies what each model really is.** The conceptual translation of every model is searched with
   `blastp` against the **complete human RefSeq proteome**, and a model counts as AGTR2 only if its best
   human match is AGTR2 itself and no other gene ties or beats it. Models that match something else are
   reported under the gene they actually match. This matters because an AGTR2 query cross-hits a long list
   of unrelated GPCRs at 27–38% identity — AGTR1, BDKRB1/2, CXCR3, CCR1/2/3/5/8, LPAR6, CMKLR2, APLNR —
   and without this step any of them can be mistaken for a faint AGTR2.
4. **Calls each species.**

   | Call | Meaning |
   |---|---|
   | `INTACT` | an AGTR2 orthologue recovering ≥ 60% of the bait in one piece, with no internal stop codon |
   | `PSEUDOGENE` | an AGTR2 orthologue carrying internal stop codons and/or recovered in several pieces |
   | `NO_AGTR2` | no AGTR2 orthologue anywhere in the assembly; the best hit is reported under the gene it belongs to |

   A stop codon counts only if it lies at least 10 residues inside its alignment, so a stop at an exon
   junction is not read as a lesion.
5. **Measures the paralogue floor.** For every assembly the highest-scoring **non**-AGTR2 model is reported
   with its gene and bit score. The dashed line in the figure is the highest of those across the run, so
   the threshold separating "AGTR2" from "a paralogue that resembles AGTR2" is **measured from the data**
   rather than chosen by eye.

## Controls, printed before the results

| Control | Requirement |
|---|---|
| **PRESENT** | every species whose AGTR2 is reported in the literature must be called `INTACT` |
| **ABSENT** | every species with a published loss must be called `NO_AGTR2` |
| **AGTR1 FOUND** | AGTR1 exists in every vertebrate, so it must be found in **every** assembly. A genome where the scan cannot find AGTR1 is a genome where the scan is not working, and its AGTR2 result means nothing. This is a control on the search itself, not on the biology. |
| **Paralogue floor** | measured, not assumed — the highest non-AGTR2 bit score in the run |
| **Search depth** | a `NO_AGTR2` call means "not among the models examined", so how far down the screen looked is reported per assembly. It passes only if the model list ran out before the `--max-models` cap (nothing was cut off) or the lowest model examined scores **below** the paralogue floor — otherwise the search stopped while a hit could still have been called AGTR2 |

The `expected_AGTR2` column of `genomes.tsv` is what defines the first two. Species marked `test` are not
used as controls.

## Scope

This analysis reports **presence or absence of an AGTR2 orthologue anywhere in the assembly**. It does not
say where the gene sits relative to its neighbours (`../TABLE1/`) nor what is left in the intergenic
interval (`../FIGURE3_Intergenic_AGTR2_PSEUDOGENE_SCreen/`).

## Running it

1. Download the human reference from NCBI Datasets, assembly **GCF_000001405.40**, choosing
   **Protein (FASTA)** and **Annotation features (GFF)**.
2. Download each assembly in `genomes.tsv` (accession column) as **Genome sequences (FASTA)** and put the
   path to each zip in the last column. `FIGURE5.py` reads the accession out of the genomic `.fna` inside
   the zip and stops if it does not match the accession listed in `genomes.tsv`.
3. Edit the paths at the top of `run_FIGURE5.sh`, then:

```bash
bash run_FIGURE5.sh
```

or directly:

```bash
python FIGURE5.py \
  --reference /path/to/human_GCF_000001405.40.zip \
  --genomes   genomes.tsv \
  --out       results \
  --db-dir    /path/to/blast_databases \
  --baits     my_baits.faa \
  --species   Homo_sapiens,Myotis_myotis,Rhinolophus_sinicus,Calypte_anna
```

`--no-figure` writes the tables and the controls and draws nothing. `run_FIGURE5.sh` ships with it switched
on: **read the data first, then draw.** Blank `NO_FIGURE` at the top of the script to produce the figure.
Re-running after the first pass is fast — the tblastn output is cached per bait in `results/raw_tblastn/`.

`--species` sets both which species are scanned and the order they appear in the figure. Without `--baits`,
the bait is extracted from the assembly named by `--bait-species` by the TABLE1 method. Only the **first**
record of the bait file is used; its length and md5 checksum appear in the generated Methods.

## Files

| File | Contents |
|---|---|
| `FIGURE5.py` | the analysis |
| `TABLE1.py` | the gene-finding library it uses (tblastn, chaining, reciprocal-best-hit orthology) |
| `genomes.tsv` | the genomes scanned, with NCBI assembly accessions and the control expectations |
| `my_baits.faa` | the AGTR2 bait; the first record is the one used |
| `run_FIGURE5.sh` | the published run, paths at the top |
| `test_FIGURE5.py` | self-test: three synthetic assemblies with a known right answer (see below) |
| `MATERIALS_AND_METHODS.py` | assembles the Methods text from the generated run files |
| `make_github_package.py` | builds the upload folder, with local paths removed |
| `results/FIGURE5.tsv` | per species: the call, and for both the best AGTR2 model and the best non-AGTR2 model — bit score, identity, coverage, E-value, number of pieces, internal stops, gene and location |
| `results/FIGURE5_all_models.tsv` | **every** chained model in every genome with the human gene it matches, so any call can be audited |
| `results/FIGURE5_controls.txt` | every control, with the species that pass or fail it |
| `results/FIGURE5_methods.txt` | Methods text generated from this run's own parameters |
| `results/FIGURE5_run_log.txt` | software versions, parameters, script checksums, assembly provenance |
| `results/FIGURE5_orthology_scan.png` / `.pdf` | the figure |

## Materials and Methods

```bash
python MATERIALS_AND_METHODS.py
```

Writes `MATERIALS_AND_METHODS.txt` from the Methods text each analysis generated about its own run, so no
sentence in it is typed by hand and none can drift from what was executed. It always includes this
analysis; if the companion analyses `TABLE1/` and `FIGURE3_Intergenic_AGTR2_PSEUDOGENE_SCreen/` have been
run it includes those too, in paper order. They are **searched for** rather than assumed to sit beside this
folder, so this directory can be moved anywhere without breaking; `--companions DIR` names the folder that
holds them if the search misses, and `$AGTR2_ROOT` does the same. `--only-this-analysis` writes the
orthology section alone, and `--out` sets the destination. Any section it cannot find is named at the end
of the run, never dropped quietly.

## Putting it on GitHub

```bash
python make_github_package.py
```

Builds `~/AGTR2_github/Figure5_AGTR2_Genomic_Orthology` — the folder to upload — from the code, the genome
list, the bait and the published results. It removes what is local to one machine: home-directory paths
become `<PATH>` in every text file (checksums, software versions and assembly provenance are kept), the
`REFERENCE` and `DB_DIR` lines of the run script become documented placeholders, and the last column of
`genomes.tsv` becomes `genomes/<accession>.zip` while every NCBI accession is kept. It writes a
`.gitignore` for the BLAST databases, the unpacked reference proteome and the tblastn cache, then prints
every file with its size and checks that your username does not appear anywhere. The output goes outside
OneDrive by default, because syncing locks files while they are being written.

Run the analysis before packaging if you want the published output included; files it cannot find are
listed rather than skipped silently.

## Self-test

```bash
python test_FIGURE5.py
```

Needs only BLAST+ and `my_baits.faa`; no genomes and no internet. It plants three synthetic assemblies
whose right answer is known and requires the program to produce it:

| Assembly | What is planted | Required call |
|---|---|---|
| `mockA` | a paralogue, plus an intact AGTR2 coding sequence, contiguous, no stops | `INTACT` |
| `mockB` | the same paralogue, plus an AGTR2 coding sequence 35% diverged, in four pieces on the minus strand 5 kb apart, each carrying an in-frame stop codon | `PSEUDOGENE` |
| `mockC` | the paralogue only | `NO_AGTR2` |

`mockC` is the test that the orthology step works: the planted paralogue is diverged enough to be a
different gene but close enough that the AGTR2 bait hits it, and the program must identify it as the
paralogue and refuse to call it AGTR2. `mockB` is the test that chaining works: each planted piece scores
about 103 bits on its own, **less** than the 202-bit contiguous paralogue alignment in the same genome, so
a best-single-HSP method would report the paralogue and call AGTR2 absent; chaining recovers the pseudogene
at 308 bits. The test also checks that the printed paralogue floor really is the highest non-AGTR2 model in
the run, that the generated Methods names the bait file and its checksum, and that no control fails.

This is a test of the program's logic, not of the biology — the sequences are synthetic.

## Notes

- BLAST databases must not sit in a cloud-synced folder (OneDrive, Dropbox) — BLAST memory-maps them and
  placeholder files cannot be mapped (`Database memory map file error`).
- Under WSL, keep `--proteome-db-dir` on the native Linux filesystem (`~/blastdb`); the LMDB component is
  unreliable on `/mnt/c` (`MDB_INVALID: File is not an LMDB file`).
- Databases built by Windows BLAST+ and by Linux BLAST+ are not interchangeable. Run the whole analysis on
  one platform.
- The `results/raw_tblastn/` cache is keyed by the bait's checksum, so changing the bait forces a fresh
  search rather than silently reusing the old one.

## Licence

MIT — see [LICENSE](LICENSE). The code is free to use, modify and redistribute, with attribution and
without warranty.

The genome assemblies analysed here are not redistributed in this repository. Every one is public and is
identified by its NCBI accession in `genomes.tsv`; download them from NCBI under their own terms.

## Citation

If you use this code, please cite the paper:

> \<AUTHORS\>. \<TITLE\>. \<JOURNAL\> (\<YEAR\>). \<DOI\>

and, if you cite the code itself rather than the paper, the archived release:

> \<ZENODO DOI\>

