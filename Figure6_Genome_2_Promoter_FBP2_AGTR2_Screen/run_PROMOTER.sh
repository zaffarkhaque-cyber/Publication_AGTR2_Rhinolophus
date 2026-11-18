#!/bin/bash
# Promoter retention screen - AGTR2 and FBP2, measured against a noise floor.
# Edit the paths below for your machine, then: bash run_PROMOTER.sh
set -e
cd "$(dirname "$0")"

REFERENCE="/path/to/human_GCF_000001405.40.zip"   # NCBI Datasets: Protein (FASTA) + Annotation features (GFF)
DB_DIR="/path/to/blast_databases"                 # NOT inside OneDrive/Dropbox
PROTEOME_DB_DIR="$HOME/blastdb"                   # under WSL keep this on the Linux filesystem
GENOMES="genomes.tsv"
QUERIES="queries.tsv"
THREADS=8

PY=python
command -v "$PY" >/dev/null 2>&1 || PY=python3

"$PY" PROMOTER.py \
  --reference        "$REFERENCE" \
  --genomes          "$GENOMES" \
  --queries          "$QUERIES" \
  --out              results \
  --db-dir           "$DB_DIR" \
  --proteome-db-dir  "$PROTEOME_DB_DIR" \
  --threads          "$THREADS"
