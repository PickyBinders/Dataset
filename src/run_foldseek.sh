#!/bin/bash
PDB_DIR=$1
OUTPUT_DIR=$2
PDB_LIST=$3
OUTPUT_NAME=$4

PDB_FS_DIR="$OUTPUT_DIR/foldseek_pdb"
if [ ! -f "$OUTPUT_DIR/foldseek_pdb"]; then
    foldseek createdb "$PDB_DIR" "$PDB_FS_DIR" --chain-name-mode 1
    foldseek createindex "$PDB_FS_DIR" tmp
fi

POCKET_FS_DIR="$OUTPUT_DIR/$OUTPUT_NAME"

awk 'NR == FNR {f[$1] = $1; next} $2 in f {print $1}' "$POCKET_FS_DIR/$OUTPUT_NAME".txt "$PDB_FS_DIR".lookup > "$POCKET_FS_DIR/$OUTPUT_NAME".tsv
foldseek createsubdb "$POCKET_FS_DIR/$OUTPUT_NAME".tsv "$PDB_FS_DIR" "$POCKET_FS_DIR/$OUTPUT_NAME" --subdb-mode 1
foldseek createsubdb "$POCKET_FS_DIR/$OUTPUT_NAME".tsv "$PDB_FS_DIR"_ss "$POCKET_FS_DIR/$OUTPUT_NAME"_ss --subdb-mode 1
foldseek createsubdb "$POCKET_FS_DIR/$OUTPUT_NAME".tsv "$PDB_FS_DIR"_ca "$POCKET_FS_DIR/$OUTPUT_NAME"_ca --subdb-mode 1

for subset_file in $POCKET_FS_DIR/subset_files/*.txt; do
    subset_name=$(basename "$subset_file" .txt)
    awk 'NR == FNR {f[$1] = $1; next} $2 in f {print $1}' $subset_file "$PDB_FS_DIR".lookup > "$POCKET_FS_DIR/subset_files/$subset_name".tsv
    foldseek createsubdb "$POCKET_FS_DIR/subset_files/$subset_name".tsv "$PDB_FS_DIR" "$POCKET_FS_DIR/subset_files/$subset_name" --subdb-mode 1
    foldseek createsubdb "$POCKET_FS_DIR/subset_files/$subset_name".tsv "$PDB_FS_DIR"_ss "$POCKET_FS_DIR/subset_files/$subset_name"_ss --subdb-mode 1
    foldseek createsubdb "$POCKET_FS_DIR/subset_files/$subset_name".tsv "$PDB_FS_DIR"_ca "$POCKET_FS_DIR/subset_files/$subset_name"_ca --subdb-mode 1

    POCKET_FS_SEARCH_DIR_SUBSET="$POCKET_FS_DIR/subset_files/search_$subset_name"

    foldseek search "$POCKET_FS_DIR/subset_files/$subset_name" "$POCKET_FS_DIR/$OUTPUT_NAME" "$POCKET_FS_SEARCH_DIR_SUBSET" tmp -a -e 0.01 --max-seqs 5000 --sort-by-structure-bits 0
    foldseek convertalis "$POCKET_FS_DIR/subset_files/$subset_name" "$POCKET_FS_DIR/$OUTPUT_NAME" "$POCKET_FS_SEARCH_DIR_SUBSET" "$POCKET_FS_DIR/subset_files/aln_$subset_name.tsv" --format-mode 4 --format-output "query,target,qlen,lddt,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln,lddtfull"
done
