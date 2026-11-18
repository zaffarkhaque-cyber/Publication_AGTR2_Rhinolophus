#!/bin/bash
# Reproduce TABLE1: AGTR2 locus evaluation across the genomes listed in genomes.tsv.
# Edit the three paths below for your machine, then: bash run_TABLE1.sh
set -e
cd "$(dirname "$0")"

REFERENCE="/path/to/human_GCF_000001405.40.zip"   # NCBI Datasets: Protein (FASTA) + Annotation features (GFF)
DB_DIR="/path/to/blast_databases"                 # NOT inside OneDrive/Dropbox
# The human proteome database is built here. Under WSL keep it on the native Linux filesystem:
# BLAST memory-maps its databases and the LMDB (.pdb) component is unreliable on a /mnt/c mount.
PROTEOME_DB_DIR="$HOME/blastdb"
THREADS=8

PY=python
command -v "$PY" >/dev/null 2>&1 || PY=python3

"$PY" TABLE1.py \
  --reference       "$REFERENCE" \
  --genomes         genomes.tsv \
  --out             results \
  --db-dir          "$DB_DIR" \
  --proteome-db-dir "$PROTEOME_DB_DIR" \
  --threads         "$THREADS"

"$PY" TABLE1_stop_codon_map.py results
"$PY" TABLE1_to_excel.py results TABLE1_results.xlsx
