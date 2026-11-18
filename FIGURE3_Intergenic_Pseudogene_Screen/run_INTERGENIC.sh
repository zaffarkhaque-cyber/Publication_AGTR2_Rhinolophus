#!/bin/bash
# Figure 3 - scan of the PLS3-SLC6A14 intergenic region for residual AGTR2 sequence.
# Edit the paths below for your machine, then: bash run_INTERGENIC.sh
set -e
cd "$(dirname "$0")"

REFERENCE="/path/to/human_GCF_000001405.40.zip"   # NCBI Datasets: Protein (FASTA) + Annotation features (GFF)
DB_DIR="/path/to/blast_databases"                 # NOT inside OneDrive/Dropbox
PROTEOME_DB_DIR="$HOME/blastdb"                   # under WSL keep this on the Linux filesystem
GENOMES="genomes.tsv"                             # in this folder; edit the last column to your downloaded zips
THREADS=8

# Species to scan, in the order they should appear in the figure.
SPECIES="Myotis_myotis,Hipposideros_larvatus,Aselliscus_stoliczkanus,Rhinolophus_ferrumequinum,Rhinolophus_sinicus,Rhinolophus_hipposideros,Rhinolophus_affinis,Apus_apus,Calypte_anna"

PY=python
command -v "$PY" >/dev/null 2>&1 || PY=python3

"$PY" INTERGENIC.py \
  --reference       "$REFERENCE" \
  --genomes         "$GENOMES" \
  --out             results \
  --db-dir          "$DB_DIR" \
  --proteome-db-dir "$PROTEOME_DB_DIR" \
  --baits           my_baits.faa \
  --species         "$SPECIES" \
  --threads         "$THREADS"
