#!/bin/bash
PDB_DIR=$1
OUTPUT_DIR=$2
AFDB_DIR=$4
OUTPUT_NAME=$5

PDB_FS_DIR="$OUTPUT_DIR/foldseek_pdb"
POCKET_FS_DIR="$OUTPUT_DIR/$OUTPUT_NAME"

AFDB_DB="$OUTPUT_DIR/afdb"
mkdir -p "$AFDB_DB"

foldseek createdb "$AFDB_DIR" "$AFDB_DB/afdb"

SEARCH_AFDB="$POCKET_FS_DIR/search_afdb"
foldseek search "$POCKET_FS_DIR/$OUTPUT_NAME" "$AFDB_DB/afdb" "$SEARCH_AFDB" tmp -a --sort-by-structure-bits 0 --min-seq-id 0.95 -c 0.9
foldseek convertalis "$POCKET_FS_DIR/$OUTPUT_NAME" "$AFDB_DB/afdb" "$SEARCH_AFDB" "$POCKET_FS_DIR/aln_afdb.tsv" --format-mode 4 --format-output "query,target,qlen,lddt,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln,lddtfull"