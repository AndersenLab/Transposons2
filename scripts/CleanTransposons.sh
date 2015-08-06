#!/bin/bash
# NOTE: shoul always run ProcessTransposons immediately prior to this script

#directories
scripts_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts
data_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data

#scripts
rename_final_results=rename_final_results.sh
resolve_contradictory=resolve_contractictory_calls.py
compatibility_script=temp_telocate_compatibility.py
recount=run_recount.sh
merge=merge_te_data.sh
reformat=reformat_te_traits.py
column_nmaes=te_column_names.py
transpose=transpose_matrix.sh
assign_cut_copy=assign_copy_or_cut_paste.py

#files
samples=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/full_sample_list.txt
results_file=/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/FINAL_RESULTS.txt
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 
repbase_fasta=/lscr2/andersenlab/kml436/repbase.fasta

# Rename results files so that they will not be overwritten in a later step
bash ${scripts_dir}/${rename_final_results} $samples
# Run compatibilty script
python ${scripts_dir}/${compatibility_script} $data_dir $samples
# Resolve contradictory calls
bash ${scripts_dir}/${resolve_contradictory_calls} $data_dir $samples
# Rename previous results files so that they are saved jsut in case but not used 
find ../results -name "FINAL*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find ../results -name "Full*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
# Recount the total totals in each strain
python ${scripts_dir}/${recount} $samples $results_file $consensus_renamed $data_dir
# Merge positions to generate the "all_nonredundant file"
bash ${scripts_dir}/${merge}
# Reformat results file
cd ../results
python ${scripts_dir}/${reformat} FINAL_RESULTS_LF.txt
# Pull correct order of column names
python ${scripts_dir}/${column_names} FINAL_RESULTS_LF.txt
# Remove trailing tab and replave it with a newline
sed  's/\t$/\n/' column_names.txt> tmp && mv tmp column_names.txt
# Add on appropriately ordered column names
cat column_names.txt FINAL_RESULTS_LF_ReF.txt >Full_Results.txt
# Transpose the matrix and add in coverage trait, remove QX strains if still present
cat Full_Results.txt | awk '$1 != “QX2265” && $1 != "QX2266" {print $0}' > tmp && mv tmp Full_Results.txt
bash ${scripts_dir}/${transpose_matrix} Full_Results.txt T_Full_Results.txt
cat T_coverage.txt |sed '1d' >tmp
cat T_Full_Results.txt tmp > all && mv all T_Full_Results.txt
# Add information on whether the TE family is a dnatransposon, retrotransposon, or unknown
python ${scripts_dir}/${assign_cut_copy} $repbase_fasta $consensus_renamed ${data_dir}/all_nonredundant.txt

