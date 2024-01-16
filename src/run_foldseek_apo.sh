#!/bin/bash
PDB_DIR=$1
OUTPUT_DIR=$2
APO_LIST=$3
OUTPUT_NAME=$4

PDB_FS_DIR="$OUTPUT_DIR/foldseek_pdb"
POCKET_FS_DIR="$OUTPUT_DIR/$OUTPUT_NAME"

APO_DB="$OUTPUT_DIR/apo"
mkdir -p "$APO_DB"

awk 'NR == FNR {f[$1] = $1; next} $2 in f {print $1}' "$APO_LIST" "$PDB_FS_DIR".lookup > "$APO_DB/apo_chains".tsv
foldseek createsubdb "$APO_DB/apo_chains".tsv "$PDB_FS_DIR" "$APO_DB/apo_chains" --subdb-mode 1
foldseek createsubdb "$APO_DB/apo_chains".tsv "$PDB_FS_DIR"_ss "$APO_DB/apo_chains_ss" --subdb-mode 1
foldseek createsubdb "$APO_DB/apo_chains".tsv "$PDB_FS_DIR"_ca "$APO_DB/apo_chains_ca" --subdb-mode 1

SEARCH_APO="$POCKET_FS_DIR/search_apo"
foldseek search "$POCKET_FS_DIR/$OUTPUT_NAME" "$APO_DIR/apo_chains" "$SEARCH_APO" tmp -a --sort-by-structure-bits 0 --min-seq-id 0.95 -c 0.9
foldseek convertalis "$POCKET_FS_DIR/$OUTPUT_NAME" "$APO_DIR/apo_chains" "$SEARCH_APO" "$POCKET_FS_DIR/aln_apo.tsv" --format-mode 4 --format-output "query,target,qlen,lddt,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln,lddtfull"