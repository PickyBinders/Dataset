#!/bin/bash

# Base directory and other common variables
FOLDSEEK_DIR=$1
LOG_DIR=$2
OUTPUT_NAME="filtered_pockets"
PDB_DIR="$FOLDSEEK_DIR/pdb/pdb"
POCKET_DIR="$FOLDSEEK_DIR/$OUTPUT_NAME"

# Main job submission
main_job_id=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH --job-name=foldseek
#SBATCH --output=$LOG_DIR/foldseek.oe
#SBATCH --mem=128G
#SBATCH --qos=30min
#SBATCH -n 64
#SBATCH -N 1

foldseek createdb /scicore/data/managed/PDB/latest/data/structures/all/mmCIF/ "$PDB_DIR" --chain-name-mode 1
foldseek createindex "$PDB_DIR" tmp_foldseek
awk 'NR == FNR {f[\$1] = \$1; next} \$2 in f {print \$1}' "$POCKET_DIR/$OUTPUT_NAME".txt "$PDB_DIR".lookup > "$POCKET_DIR/$OUTPUT_NAME".tsv
foldseek createsubdb "$POCKET_DIR/$OUTPUT_NAME".tsv "$PDB_DIR" "$POCKET_DIR/$OUTPUT_NAME" --subdb-mode 1
foldseek createsubdb "$POCKET_DIR/$OUTPUT_NAME".tsv "$PDB_DIR"_ss "$POCKET_DIR/$OUTPUT_NAME"_ss --subdb-mode 1
foldseek createsubdb "$POCKET_DIR/$OUTPUT_NAME".tsv "$PDB_DIR"_ca "$POCKET_DIR/$OUTPUT_NAME"_ca --subdb-mode 1
EOF
)

if [ -z "$main_job_id" ]; then
    echo "Failed to submit the main job or capture its ID."
    exit 1
fi

sbatch --dependency=afterok:$main_job_id << EOF
#!/bin/bash
#SBATCH --job-name=foldseek_subsets
#SBATCH --output="$LOG_DIR"/foldseek_subsets.oe
#SBATCH --mem=256G
#SBATCH --qos=6hours
#SBATCH -n 128

for subset_file in "$POCKET_DIR/subset_files/"*.txt; do
    subset_name=\$(basename "\$subset_file" .txt)
    awk 'NR == FNR {f[\$1] = \$1; next} \$2 in f {print \$1}' \$subset_file "$PDB_DIR".lookup > "$POCKET_DIR/subset_files/\$subset_name".tsv
    foldseek createsubdb "$POCKET_DIR/subset_files/\$subset_name".tsv "$PDB_DIR" "$POCKET_DIR/subset_files/\$subset_name" --subdb-mode 1
    foldseek createsubdb "$POCKET_DIR/subset_files/\$subset_name".tsv "$PDB_DIR"_ss "$POCKET_DIR/subset_files/\$subset_name"_ss --subdb-mode 1
    foldseek createsubdb "$POCKET_DIR/subset_files/\$subset_name".tsv "$PDB_DIR"_ca "$POCKET_DIR/subset_files/\$subset_name"_ca --subdb-mode 1

    POCKET_SEARCH_DIR_SUBSET="$POCKET_DIR/subset_files/search_\$subset_name"

    foldseek search "$POCKET_DIR/subset_files/\$subset_name" "$POCKET_DIR/$OUTPUT_NAME" "\$POCKET_SEARCH_DIR_SUBSET" tmp_foldseek -a -e 0.01 --max-seqs 5000 --sort-by-structure-bits 0
    foldseek convertalis "$POCKET_DIR/subset_files/\$subset_name" "$POCKET_DIR/$OUTPUT_NAME" "\$POCKET_SEARCH_DIR_SUBSET" "$POCKET_DIR/subset_files/aln_\$subset_name.tsv" --format-mode 4 --format-output "query,target,qlen,lddt,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln,lddtfull"
done
EOF
