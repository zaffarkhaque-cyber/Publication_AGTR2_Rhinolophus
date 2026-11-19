# Materials and Methods

Generated from the run itself by `code/make_package.py`. Every number below is read from the
manifests and the run log, not transcribed.

## Genomes

23 genome assemblies were searched (16 bat species and one outgroup, *Homo sapiens*, plus bird
genomes used only for clade-matched reference assignment). Each was taken from an NCBI datasets
archive; the accession recorded for each in `manifests/genomes.tsv` was read from inside the archive
rather than typed, so it cannot drift from the file that was searched.

## Queries

13 protein queries and 18 promoter queries were taken from a single hand-written design file
(`manifests/DESIGN_PPEF2_NOX5.txt`) and parsed by `code/parse_design.py`, which writes the manifests
and one FASTA per query. The md5 of the design file is recorded in `manifests/genomes.tsv`, and the
md5 of every query is recorded in `results/RUN_LOG.txt`. No sequence was retyped at any point: each
query is a copy of a block in the design file, and three promoter queries are a declared suffix of
their block (the last 600 bases, the transcription start site being the 3' end), recorded in
`manifests/promoters.tsv` alongside the untrimmed block length and md5.

## Gene search

Genes were sought with `tblastn` and the resulting high-scoring pairs chained per scaffold and
strand by dynamic programming, with a maximum intron of 1,000,000 bp, so that a gene broken into
pieces is recovered as one model rather than lost. Every candidate chain was examined, not only the
best-scoring one, because a paralogue routinely outscores the gene: AGTR1 over AGTR2, PPEF1 over
PPEF2, IL1RN over IL1F10.

Each model was identified by reciprocal best hit against complete reference proteomes, and a model
counts as the query gene only when its best match is that gene and nothing ties or beats it. Models
were judged against a reference from their own clade, because a mammal proteome names some bird
genes after their paralogues. NCBI readthrough annotations (hyphenated symbols such as SPESP1-NOX5)
are treated as the gene rather than as a rival.

Calls: INTACT at or above 0.60, min stops 2 query coverage with fewer than 2 confident stop codons;
PSEUDOGENE at or above that coverage with 2 or more; FRAGMENT below it; AMBIGUOUS where
same-clade references disagree; ABSENT where no chain is confirmed as the gene. A stop codon counts
as confident only when it lies at least 10 residues inside an exon hit and has at least 50% identity
over the 15 residues on each side, which excludes stops produced by an alignment running into an
intron.

## Assembly integrity at each locus

For every call, the genomic interval of the model plus 20,000 bp on each side was extracted from the
BLAST database and its ambiguous bases counted (`code/locus_gap_check.py`, output in
`results/LOCUS_GAPS.tsv`). An absence means something only where that count is zero: continuous
sequence was assembled through the region and the gene is not in it.

The following genome(s) carry ambiguous sequence at these loci and are therefore not assessable; they are searched and reported in the tables but excluded from the figures and from every claim: *Hipposideros armiger*, *Pteropus alecto*.

## Controls

Three controls run on every analysis and are reported in `results/CONTROLS.txt`:

- **Self-hit.** Every query must be recovered from the genome it was taken from. All protein queries
  returned INTACT and all promoter queries returned DETECTED in their own source genome.
- **Control genes.** SELE and ACE carry no claim in this work. A control gene is assessed only in a
  clade it has a query from, and counts as recovered in a genome if any of its queries returns it
  intact. SELE was recovered intact in all 23 genomes searched; ACE in all 17 mammals.
- **Clade cover.** Every clade carrying a claim has a control query from that clade, so a negative
  can be distinguished from the assay failing in that clade.

## Repeat accounting

Promoter searches report, in addition to aligned length, the number of separate genomic loci each
query base is recovered from, and the fraction of the query recovered from sequence held at fewer
than 10 loci. Low-complexity filtering (DUST) remains enabled throughout. This distinction is not
cosmetic: 45 of 414 promoter cells align well above the shuffled-sequence noise floor with excellent
E-values while containing no unique sequence at all, and would be false positives on aligned length
alone.

## Reproducibility

`code/test_interrogate.py` contains the control tests for the functions that decide a call, each a
positive case paired with a negative. `MANIFEST.md5` lists a checksum for every file in this
package. `command: interrogate.py --mammal-reference <PATH>/OneDrive/Desktop/BATCH GENOMES/Group 12/COMPLETE_Homo-Sapiens_ncbi_dataset.zip --mammal-reference2 <PATH>/Downloads/Rousettus_aegyptiacius_whole_ncbi_dataset.zip --bird-reference <PATH>/Downloads/WHOLE_Apus_apus_ncbi_dataset.zip --db-dir <PATH>/bat_promoter_blast_databases --out <PATH>/OneDrive/Desktop/NOX5_PPEF2_Genome_Interrogation/results`

promoter call needs bits > floor, E < 1e-05, and >= 10% of the query recovered from sequence held at fewer than 10 loci
