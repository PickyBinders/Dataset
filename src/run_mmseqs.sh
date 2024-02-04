#!/bin/bash

# Base directory and other common variables
MMSEQS_DIR=$1
LOG_DIR=$2
OUTPUT_NAME="filtered_pockets"
PDB_DIR="$MMSEQS_DIR/pdb/pdb"
POCKET_DIR="$MMSEQS_DIR/$OUTPUT_NAME"

# Main job submission
main_job_id=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH --job-name=mmseqs
#SBATCH --output=$LOG_DIR/mmseqs.oe
#SBATCH --mem=128G
#SBATCH --qos=30min
#SBATCH -n 64
#SBATCH -N 1

mmseqs createdb /scicore/data/managed/PDB/latest/derived_data/pdb_seqres.txt "$PDB_DIR"
mmseqs createindex "$PDB_DIR" tmp_mmseqs
awk 'NR == FNR {f[\$1] = \$1; next} \$2 in f {print \$1}' "$POCKET_DIR/$OUTPUT_NAME".txt "$PDB_DIR".lookup > "$POCKET_DIR/$OUTPUT_NAME".tsv
mmseqs createsubdb "$POCKET_DIR/$OUTPUT_NAME".tsv "$PDB_DIR" "$POCKET_DIR/$OUTPUT_NAME" 
EOF
)

if [ -z "$main_job_id" ]; then
    echo "Failed to submit the main job or capture its ID."
    exit 1
fi

sbatch --dependency=afterok:$main_job_id << EOF
#!/bin/bash
#SBATCH --job-name=mmseqs_subsets
#SBATCH --output="$LOG_DIR"/mmseqs_subsets.oe
#SBATCH --mem=256G
#SBATCH --qos=6hours
#SBATCH -n 64

POCKET_SEARCH_DIR="$POCKET_DIR/search"
mmseqs search "$POCKET_DIR/$OUTPUT_NAME" "$POCKET_DIR/$OUTPUT_NAME" "\$POCKET_SEARCH_DIR" tmp_mmseqs -a -e 0.01 --max-seqs 5000
mmseqs convertalis "$POCKET_DIR/$OUTPUT_NAME" "$POCKET_DIR/$OUTPUT_NAME" "\$POCKET_SEARCH_DIR" "$POCKET_DIR/aln.tsv" --format-mode 4 --format-output "query,target,qlen,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qaln,taln"
EOF
