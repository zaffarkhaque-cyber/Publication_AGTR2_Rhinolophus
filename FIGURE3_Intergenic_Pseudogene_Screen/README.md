# Figure 3 — intergenic AGTR2 pseudogene screen

**Question:** in a genome where AGTR2 is not found, is there a decayed remnant (pseudogene) left
between the flanking genes PLS3 and SLC6A14, or is the locus deleted outright?

**This directory is self-contained.** Everything needed to reproduce the figure is here: the
analysis, the gene-finding library it uses, the genome list, the bait sequences, the run script and
the output of the published run.

The companion directory `TABLE1` answers the separate question of whether AGTR2 is present anywhere
in each genome. `TABLE1.py` is shared between the two analyses, and the copy in the `TABLE1`
directory is the master — if you change it there, copy it here as well.

## What the analysis does

1. **Locates the interval.** The human PLS3 and SLC6A14 proteins are searched against each assembly
   with `tblastn`, hits are chained into collinear gene models, and orthology is assigned by
   reciprocal best hit against the complete human proteome. Where both orthologues lie on one
   scaffold, the DNA between them is read directly from the assembly — no browser downloads, so
   every interval carries an accession and coordinates.
2. **Searches it.** Each AGTR2 protein bait is searched against each interval with `tblastn`
   (protein against translated DNA), which stays sensitive to diverged remnants long after
   nucleotide comparison fails.
3. **Measures the noise floor.** Each interval is shuffled 20 times, preserving length and base
   composition, and rescanned. The floor is the greater of the highest shuffled score and the mean
   plus three standard deviations.
4. **Measures the detection limit.** Each bait's coding sequence is mutated at 10–60%, embedded in
   shuffled sequence and put through the same scan — so a negative result states what it excludes.

A hit counts as a remnant only if it exceeds the floor **and** has E < 0.001. Either test alone
admits chance alignments in a long interval.

## Controls, printed before the results

| Control | Requirement |
|---|---|
| **INTACT** | species whose AGTR2 lies in the interval intact must be found |
| **REMNANT** | species whose AGTR2 in the interval is a **pseudogene** must still be found — without this a negative result is uninterpretable |
| **Noise floor** | measured from shuffled sequence, per species, per bait |
| **Detection limit** | measured, not assumed |
| **Assembly gaps** | N content and longest N run per interval, so an absence inside a gap is not read as a deletion |

## Scope

This analysis reports **only** what lies between PLS3 and SLC6A14. A species may carry AGTR2
sequence elsewhere in its genome; `../TABLE1/` is what answers that.

## Running it

1. Download the human reference from NCBI Datasets, assembly **GCF_000001405.40**, choosing
   **Protein (FASTA)** and **Annotation features (GFF)**.
2. Download each assembly in `genomes.tsv` (accession column) as **Genome sequences (FASTA)** and
   put the path to each zip in the last column.
3. Edit the paths at the top of `run_INTERGENIC.sh`, then:

```bash
bash run_INTERGENIC.sh
```

or directly:

```bash
python INTERGENIC.py \
  --reference /path/to/human_GCF_000001405.40.zip \
  --genomes   genomes.tsv \
  --out       results \
  --db-dir    /path/to/blast_databases \
  --baits     my_baits.faa \
  --species   Myotis_myotis,Rhinolophus_sinicus,Apus_apus
```

`--species` sets both which species are scanned and the order they appear in the figure.
Without `--baits`, baits are taken from the assemblies named by `--bait-species`.

## Files

| File | Contents |
|---|---|
| `INTERGENIC.py` | the analysis |
| `TABLE1.py` | the gene-finding library it uses (tblastn, chaining, reciprocal-best-hit orthology) |
| `genomes.tsv` | the genomes scanned, with NCBI assembly accessions |
| `run_INTERGENIC.sh` | the published run, paths at the top |
| `my_baits.faa` | the AGTR2 bait proteins used; their lengths and md5 checksums appear in the generated Methods |
| `results/INTERGENIC.tsv` | per species: interval coordinates, length, N content, and per bait the identical residues, bit score, E-value, noise floor and shuffled statistics |
| `results/INTERGENIC_detection_limit.tsv` | score against bait sequences mutated at 10–60% |
| `results/INTERGENIC_controls.txt` | every control, with the species that pass or fail it |
| `results/INTERGENIC_methods.txt` | Methods text generated from this run's own parameters |
| `results/INTERGENIC_run_log.txt` | software versions, parameters, input checksums, interval coordinates |
| `results/INTERGENIC_scan.png` / `.pdf` | the figure |

## Notes

- BLAST databases must not sit in a cloud-synced folder (OneDrive, Dropbox) — BLAST memory-maps
  them and placeholder files cannot be mapped (`Database memory map file error`).
- Under WSL, keep `--proteome-db-dir` on the native Linux filesystem (`~/blastdb`); the LMDB
  component is unreliable on `/mnt/c` (`MDB_INVALID: File is not an LMDB file`).
