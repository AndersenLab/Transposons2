#!/bin/bash
# NOTE: shoul always run ProcessTransposons immediately prior to this script
# USE: CleanTransposons.sh (run in the data dir)

#directories
scripts_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts
data_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data

#scripts
rename_final_results=rename_final_results.sh
resolve_contradictory=resolve_contradictory_calls.py
compatibility_script=temp_telocate_compatibility.py
recount=run_recount.sh
merge=merge_te_data.sh
reformat=reformat_te_traits.py
column_names=te_column_names.py
transpose=transpose_matrix.sh
assign_cut_copy=assign_copy_or_cut_paste.py
coverage=add_coverage_data.py

#files
samples=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/full_sample_list.txt
results_file=/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/FINAL_RESULTS.txt
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 
repbase_fasta=/lscr2/andersenlab/kml436/repbase.fasta
bam_stats=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/eav.txt
# Rename results files so that they will not be overwritten in a later step
bash ${scripts_dir}/${rename_final_results} $samples
# Run compatibilty script
cd ${scripts_dir}/../results/
python ${scripts_dir}/${compatibility_script} $data_dir $samples
# Resolve contradictory calls
cd $data_dir
python ${scripts_dir}/${resolve_contradictory} $data_dir $samples
# Rename previous results files so that they are saved jsut in case but not used 
find ../results -name "FINAL*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find ../results -name "Full*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
# Recount the total totals in each strain
bash ${scripts_dir}/${recount} $samples $results_file $consensus_renamed $data_dir
# Merge positions to generate the "all_nonredundant file"
bash ${scripts_dir}/${merge}
# Reformat results file
cd ../results
python ${scripts_dir}/${reformat} FINAL_RESULTS_LF.txt
# Pull correct order of column names
python ${scripts_dir}/${column_names} FINAL_RESULTS_LF.txt
# Remove trailing tab and replave it with a newline
sed  's/\t$/\n/' column_names.txt> tmp && mv tmp column_names.txt
# Add on appropriately ordered column names to FINAL_RESULTS_LF_REF
cat column_names.txt FINAL_RESULTS_LF_ReF.txt >Full_Results.txt
# Add on appropriately ordered column names to FINAL_RESULTS.txt
cat te_names_ordered.txt FINAL_RESULTS.txt > tmp && mv tmp FINAL_RESULTS.txt
# Transpose the matrix and add in coverage trait, remove QX strains if still present
cat Full_Results.txt | awk '$1 != "QX2265" && $1 != "QX2266" {print $0}' > tmp && mv tmp Full_Results.txt
bash ${scripts_dir}/${transpose} Full_Results.txt T_Full_Results.txt
python ${scripts_dir}/${coverage} Full_Results.txt $bam_stats
cat coverage_and_te_counts.txt |sed '1d'| cut -f1,5 > temp_cov
bash ../scripts/transpose_matrix.sh temp_cov temp_cov2
cat temp_cov2 | awk 'NR==1 {print "trait\t"$0}; NR==2 {print "coverage\t"$0}' > T_coverage.txt
cat T_coverage.txt |sed '1d' >tmp
cat T_Full_Results.txt tmp > all && mv all T_Full_Results.txt
# Add information on whether the TE family is a dnatransposon, retrotransposon, or unknown
python ${scripts_dir}/${assign_cut_copy} $repbase_fasta $consensus_renamed ${data_dir}/all_nonredundant.txt
# remove lines with all zeros (no instance of that transposition event in any of the samples)
mv T_Full_Results.txt original_T_Full_Results.txt
cat original_T_Full_Results.txt | awk '{for (i=2;i<=NF;i++) {if($i !=0){print $0;break}}}' > T_Full_Results.txt

