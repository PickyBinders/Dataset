validation_data_dir=$1
log_dir=$2

mkdir -p "$validation_data_dir"

num_commands=0
CMD_FILE=extract_validation_data_commands.cmd
rm -f $CMD_FILE
for folder in /scicore/data/managed/PDB/latest/validation_reports/*; do
    folder_stem=${folder##*/}
    results_file="$validation_data_dir/Results_$folder_stem.csv"
    echo "PDBValidationExtract" -v -n 10 -d 6 "$results_file" "$log_dir" "/scicore/data/managed/PDB/latest/validation_reports/$folder_stem/*/*_validation.xml.gz" >> $CMD_FILE
    num_commands=$((num_commands+1))
done

log_dir="$log_dir/array_logs"
mkdir -p $log_dir

cat << EOF | sbatch
#!/bin/bash
#SBATCH --job-name=extract_validation_data
#SBATCH --qos=1day
#SBATCH --output="$log_dir"/%a.out
#SBATCH --mem=32G
#SBATCH -n 128
#SBATCH --array=1-"$num_commands"

export PDB_CIF_PATH="/scicore/data/managed/PDB/latest/data/structures/divided/mmCIF/{short}/{full}.cif.gz"
export PDB_VALIDATION_PATH="/scicore/data/managed/PDB/latest/validation_reports/{short}/{full}/{full}_validation.xml.gz"
ml OpenStructure

\$(head -\$SLURM_ARRAY_TASK_ID $CMD_FILE | tail -1)
EOF