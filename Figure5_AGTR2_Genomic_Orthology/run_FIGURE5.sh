#!/bin/bash
# Figure 5 - whole-genome AGTR2 orthology screen.
# Edit the paths below for your machine, then: bash run_FIGURE5.sh
set -e
cd "$(dirname "$0")"

REFERENCE="/path/to/human_GCF_000001405.40.zip"   # NCBI Datasets: Protein (FASTA) + Annotation features (GFF)
DB_DIR="/path/to/blast_databases"                 # NOT inside OneDrive/Dropbox
PROTEOME_DB_DIR="$HOME/blastdb"                   # under WSL keep this on the Linux filesystem
GENOMES="genomes.tsv"                             # in this folder; edit the last column to your downloaded zips
THREADS=8

# Species to scan, in the order they should appear in the figure.
# Homo sapiens / Myotis / Pteropus / Megaderma are AGTR2-present controls (Megaderma is the rhinolophoid
# outgroup); Apus is the bird control; the three hummingbirds are the published-loss control.
SPECIES="Homo_sapiens,Myotis_myotis,Pteropus_alecto,Megaderma_spasma,Hipposideros_armiger,Hipposideros_larvatus,Aselliscus_stoliczkanus,Rhinolophus_ferrumequinum,Rhinolophus_sinicus,Rhinolophus_hipposideros,Rhinolophus_affinis,Apus_apus,Archilochus_colubris,Calypte_anna,Heliangelus_exortis"

# The data was read and the controls checked first, so the figure is now drawn.
# Put --no-figure back here to go back to tables only.
NO_FIGURE=""

# Optional: shade intact orthologues below this % identity paler, with the cutoff named in the legend.
# Presentation only - no call and no number changes. Uncomment to separate Apus apus (58.0%) from the
# 93.6-100% orthologues, as in the earlier hand-made figure.
# DIVERGENT="--divergent-below 70"
DIVERGENT=""

PY=python
command -v "$PY" >/dev/null 2>&1 || PY=python3

"$PY" FIGURE5.py $NO_FIGURE $DIVERGENT \
  --reference       "$REFERENCE" \
  --genomes         "$GENOMES" \
  --out             results \
  --db-dir          "$DB_DIR" \
  --proteome-db-dir "$PROTEOME_DB_DIR" \
  --baits           my_baits.faa \
  --species         "$SPECIES" \
  --threads         "$THREADS"
