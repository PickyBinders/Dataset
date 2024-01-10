cif_data_dir=$1

num_commands=1
rm extract_cif_data_commands.cmd
for folder in /scicore/data/managed/PDB/latest/data/structures/divided/mmCIF/*; do
    folder_stem=${folder##*/}
    results_file="$cif_data_dir/$folder_stem.json"
    echo python extract_cif_data.py $folder $results_file --ignore 7Y7A >> extract_cif_data_commands.cmd
    num_commands=$((num_commands+1))
done

log_dir=$2
mkdir -p $log_dir

# Use a heredoc to create the script
cat << EOF | sbatch
#!/bin/bash
#SBATCH --job-name=extract_cif_data
#SBATCH --qos=6hours
#SBATCH --output="$log_dir"/%a.out
#SBATCH --mem=32G
#SBATCH --array=1-"$num_commands"

ml OpenStructure
\$(head -\$SLURM_ARRAY_TASK_ID extract_cif_data_commands.cmd | tail -1)
EOF