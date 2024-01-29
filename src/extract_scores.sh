df_file=$1
cif_data_dir=$2
output_dir=$3
foldseek_aln_dir=$4
mmseqs_aln_dir=$5

num_commands=0
CMD_FILE=extract_scores_commands.cmd
rm -f $CMD_FILE
for folder in /scicore/data/managed/PDB/latest/data/structures/divided/mmCIF/*; do
    folder_stem=${folder##*/}
    echo python extract_scores.py --df_file $df_file --cif_data_dir $cif_data_dir --aln_file_foldseek "$foldseek_aln_dir/aln_$folder_stem.tsv" --aln_file_mmseqs "$mmseqs_aln_dir/aln_$folder_stem.tsv" --score_file "$output_dir/$folder_stem"_scores.tsv --overwrite >> $CMD_FILE
    num_commands=$((num_commands+1))
done

log_dir=$6
mkdir -p $log_dir

cat << EOF | sbatch
#!/bin/bash
#SBATCH --job-name=extract_scores
#SBATCH --qos=6hours
#SBATCH --output="$log_dir"/%a.out
#SBATCH --mem=32G
#SBATCH --array=1-"$num_commands"

\$(head -\$SLURM_ARRAY_TASK_ID $CMD_FILE | tail -1)
EOF