df_file=$1
cif_data_dir=$2
aln_dir=$3
score_dir=$4

num_commands=1
rm -f extract_scores_commands.cmd
for folder in /scicore/data/managed/PDB/latest/data/structures/divided/mmCIF/*; do
    folder_stem=${folder##*/}
    echo python similarity_scoring.py --df_file $df_file --cif_data_dir $cif_data_dir --prefix $folder_stem --aln_dir $aln_dir --score_dir $score_dir --overwrite >> extract_scores_commands.cmd
    num_commands=$((num_commands+1))
done

log_dir=$5
mkdir -p $log_dir

# Use a heredoc to create the script
cat << EOF | sbatch
#!/bin/bash
#SBATCH --job-name=extract_scores
#SBATCH --qos=6hours
#SBATCH --output="$log_dir"/%a.out
#SBATCH --mem=32G
#SBATCH --array=1-"$num_commands"

\$(head -\$SLURM_ARRAY_TASK_ID extract_scores_commands.cmd | tail -1)
EOF