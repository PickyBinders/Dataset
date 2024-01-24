score_dir=$1
output_dir=$2
score_names=( "protein_lddt" "protein_lddt_qcov" "protein_qcov" "protein_fident" "pocket_lddt" "pocket_lddt_qcov" "pocket_qcov" "pocket_fident" "pli_qcov" "pli_qcov_residue" "pli_qcov_nowater" )
suffixes=( "_weighted_sum" "_weighted_max" "_max" )
num_commands=1
rm -f graph_clustering_commands.cmd
for score_name in "${score_names[@]}"; do
    for suffix in "${suffixes[@]}"; do
        echo python graph_clustering.py --score_dir $score_dir --score_name $score_name$suffix --output_dir $output_dir >> graph_clustering_commands.cmd
        num_commands=$((num_commands+1))
done

log_dir=$3
mkdir -p $log_dir

cat << EOF | sbatch
#!/bin/bash
#SBATCH --job-name=graph_clustering
#SBATCH --qos=6hours
#SBATCH --output="$log_dir"/%a.out
#SBATCH --mem=128G
#SBATCH --array=1-"$num_commands"

\$(head -\$SLURM_ARRAY_TASK_ID graph_clustering_commands.cmd | tail -1)
EOF