# AGTR2 locus analyses

Two linked analyses of the AGTR2 locus across vertebrate genome assemblies.

**1. `TABLE1.py` — is AGTR2 present in the genome?**
Locates AGTR2 and its X-chromosome neighbours PLS3 and SLC6A14 (plus the paralogue AGTR1) in each
assembly, measures the PLS3–SLC6A14 interval and its N content, and where AGTR2 is present but
disrupted reports every stop codon and reading-frame change on human AGTR2 coordinates, including
those shared between species.

**2. `INTERGENIC.py` — is there any AGTR2 sequence left between PLS3 and SLC6A14?**
Reads the intergenic DNA directly from each assembly and searches it with AGTR2 protein baits, with
a noise floor measured from shuffled sequence and a measured detection limit, so that a negative
result states what it excludes.

Neither uses genome annotation; the assemblies are searched directly, so unannotated and newly
released genomes can be used.

## Why the results can be trusted

Both analyses report their controls before their results, and both fail visibly rather than quietly.

| Control | What it requires |
|---|---|
| **Gene presence** | species with AGTR2 reported in the literature must return an intact gene; species with published loss must return none |
| **Disruption calling** | every species with an intact AGTR2 must show *no* stop codon and *no* frame change — if intact genes carry lesions, lesion calls cannot be trusted |
| **Scan sensitivity** | a species whose AGTR2 is a pseudogene must still be detected inside its interval — without this a negative result is uninterpretable |
| **Noise floor** | measured by shuffling each interval, not assumed |
| **Detection limit** | measured by mutating the bait at increasing rates — states how diverged a remnant could be and still be found |
| **Assembly gaps** | N content and longest N run reported per interval, so an absence inside a gap is not mistaken for a deletion |

The orthology step is what separates AGTR2 from AGTR1 and real PLS3 from PLS1, LCP1 and retrocopies:
every gene model, and every exon hit within it, must have the query gene itself as its best match in
the complete human proteome.

## Requirements

- Python 3.8+ (standard library only for `TABLE1.py` and `INTERGENIC.py`)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/) 2.13+ (`tblastn`, `blastp`, `makeblastdb`, `blastdbcmd`) on the `PATH`
- Optional, for figures and the workbook: `matplotlib`, `pandas`, `openpyxl`

## Quick start

1. Download the human reference from NCBI Datasets, assembly **GCF_000001405.40**, selecting
   **Protein (FASTA)** and **Annotation features (GFF)**.
2. Download each assembly listed in `genomes.tsv` (accession column) as **Genome sequences (FASTA)**
   and put the path to each zip in the last column.
3. Edit the paths at the top of `run_TABLE1.sh` and `run_INTERGENIC.sh`, then:

```bash
bash run_TABLE1.sh        # Table 1, the lesion table, the stop-codon figure, the workbook
bash run_INTERGENIC.sh    # the intergenic scan and its figure
python MATERIALS_AND_METHODS.py   # assembles the Methods from what the runs recorded
```

Or call the scripts directly:

```bash
python TABLE1.py --reference human.zip --genomes genomes.tsv --out results \
                 --db-dir /path/to/blast_databases

python INTERGENIC.py --reference human.zip --genomes genomes.tsv --out results_intergenic \
                     --db-dir /path/to/blast_databases --baits my_baits.faa \
                     --species Myotis_myotis,Rhinolophus_sinicus,Apus_apus
```

BLAST databases are built once in `--db-dir` and reused, and only if built from the same file, so a
swapped assembly cannot be used silently. `INTERGENIC.py` also checks that each genome file really
contains the accession `genomes.tsv` claims for it.

## Inputs

**`genomes.tsv`** — tab-separated, `#` comments ignored:

```
species <TAB> group <TAB> expected_AGTR2 <TAB> assembly accession <TAB> path to the NCBI zip
```

`expected_AGTR2` is `present` (control: must be found), `absent` (published loss; control: must not
be found) or `test`. The accession is permanent; the path is machine-specific and must be edited.

**`my_baits.faa`** — optional protein FASTA of your own AGTR2 bait sequences for `INTERGENIC.py`.
Each header becomes a bar group in the figure, and each bait's length and md5 are recorded in the
generated Methods, so the figure is tied to the exact sequences used. Without `--baits`, baits are
taken from the assemblies of `--bait-species`.

## Outputs

`TABLE1.py` → `results/`

| File | Contents |
|---|---|
| `TABLE1.tsv` | one row per genome: AGTR2 call, locations, interval length, N content, best remaining hit |
| `TABLE1_disruptions.tsv` | every stop codon and frame change, with codon and genome coordinate |
| `TABLE1_shared_disruptions.txt` | which lesions are shared between species; disruption negative control |
| `TABLE1_gene_models.tsv` | every gene model for every query and genome, with its orthology assignment |
| `TABLE1_alignment.fasta` | each species' AGTR2 projected onto human coordinates |
| `TABLE1_controls.txt`, `TABLE1_methods.txt`, `TABLE1_run_log.txt` | controls, generated Methods, versions and checksums |

`INTERGENIC.py` → `results_intergenic/`

| File | Contents |
|---|---|
| `INTERGENIC.tsv` | per species: interval coordinates, length, N content, and per bait the identical residues, bit score, E-value, noise floor and shuffled statistics |
| `INTERGENIC_detection_limit.tsv` | score against bait sequences mutated at 10–60% |
| `INTERGENIC_scan.png` / `.pdf` | the figure |
| `INTERGENIC_controls.txt`, `INTERGENIC_methods.txt`, `INTERGENIC_run_log.txt` | controls, generated Methods, versions and checksums |

`MATERIALS_AND_METHODS.py` combines both generated Methods into `MATERIALS_AND_METHODS.txt`. Nothing
in it is typed by hand: every figure is read from the run that produced it, so the text cannot drift
from what was executed.

## Reproducing the published runs

The `results/` and `results_intergenic/` directories hold the output of the published runs. With the
same assemblies and the same commands, every file should be reproduced apart from timestamps and
paths. Each run log records the script checksum, BLAST version and all parameters.

## Notes

- **BLAST databases must not live on a cloud-synced folder** (OneDrive, Dropbox). BLAST memory-maps
  its databases and placeholder files cannot be mapped; the symptom is `Database memory map file error`.
- **Under WSL**, use `--proteome-db-dir` to build the human proteome database on the native Linux
  filesystem (for example `~/blastdb`). Its LMDB component is unreliable on a `/mnt/c` mount; the
  symptom is `MDB_INVALID: File is not an LMDB file`.
- Chromosome names come from each assembly's own FASTA titles, so numbering is not comparable between
  species. Only the X is meaningful across mammals; in birds the orthologous genes are autosomal.
- An absence is only as strong as the interval evidence. An interval containing a long run of N cannot
  exclude a hidden gene, which is why the N columns are reported per species.
- `INTERGENIC.py` reports only what lies **between** PLS3 and SLC6A14. A species may carry AGTR2
  sequence elsewhere in its genome; `TABLE1.py` is what answers that.

## Licence

TODO — add a licence before publishing (MIT is a common choice for analysis code).

## Citation

TODO — add the paper reference once available.
