score_dir=$1
output_dir=$2

FOLDSEEK_SCORE_NAMES=( "protein_lddt" "protein_lddt_qcov" "protein_qcov" "protein_fident" "protein_fident_qcov" "pocket_lddt" "pocket_lddt_qcov" "pocket_qcov" "pocket_fident" "pocket_fident_qcov" "pli_qcov" "pli_qcov_residue" "pli_qcov_nowater" )
MMSEQS_SCORE_NAMES=("protein_qcov" "protein_fident" "protein_fident_qcov" "pocket_qcov" "pocket_fident" "pocket_fident_qcov" "pli_qcov" "pli_qcov_residue" "pli_qcov_nowater" )

suffixes=( "_weighted_sum" "_weighted_max" "_max" )
types=( "strong" "weak" )
num_commands=0
CMD_FILE=graph_clustering_commands.cmd
rm -f $CMD_FILE
for score_name in "${FOLDSEEK_SCORE_NAMES[@]}"; do
    for suffix in "${suffixes[@]}"; do
        for type in "${types[@]}"; do
            echo python graph_clustering.py --score_dir $score_dir --score_name "$score_name"_foldseek"$suffix" --output_dir $output_dir --type $type >> $CMD_FILE
            num_commands=$((num_commands+1))
        done
    done
done
for score_name in "${MMSEQS_SCORE_NAMES[@]}"; do
    for suffix in "${suffixes[@]}"; do
        for type in "${types[@]}"; do
            echo python graph_clustering.py --score_dir $score_dir --score_name "$score_name"_mmseqs"$suffix" --output_dir $output_dir --type $type >> $CMD_FILE
            num_commands=$((num_commands+1))
        done
    done
done

log_dir=$3
mkdir -p $log_dir

cat << EOF | sbatch
#!/bin/bash
#SBATCH --job-name=graph_clustering
#SBATCH --qos=6hours
#SBATCH --output="$log_dir"/%a.out
#SBATCH --mem=384G
#SBATCH --array=1-"$num_commands"

\$(head -\$SLURM_ARRAY_TASK_ID $CMD_FILE | tail -1)
EOF