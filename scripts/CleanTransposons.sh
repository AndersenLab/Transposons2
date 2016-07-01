#!/bin/bash
# NOTE: should always run ProcessTransposons immediately prior to this script
# USE: CleanTransposons.sh (run in the data dir)

#directories
scripts_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts
data_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data
results_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/results

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
assign_cut_copy_R2=assign_copy_or_cut_paste_round2.py
coverage=add_coverage_data.py
fix_cov=fix_cov.py
NA=check_NAs.sh
kin=kin_hash.py
kin_AF=kin_hash_AF.py
frac_matrix_ins=te_totals_frac_ins.py
frac_matrix_AF=te_totals_frac_AF.py
NA_pos=NA_per_ins.py
reduce_ins=reduce_insertion_positions.py
frac_matrix_ins_reduced=te_totals_frac_ins_reduced.py
activity=activity_calculator.py
CtCp=generate_CtCp.py
gene_interrupt=GENE.sh
cerfinder=finds_cers.py
count_classes=total_class.py
count_classes_reduced=total_class_reduced.py
find_outliers=outliers.py
interest=pull_strains_ins_info.py
id=give_id.py
id_reduced=give_id_reduced.py
fasta_supp=fasta_supp_table.py

#files
samples=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/master_sample_list.txt
#samples=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/test_list.txt
results_file=/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/FINAL_RESULTS.txt
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 
repbase_fasta=/lscr2/andersenlab/kml436/repbase.fasta
bam_stats=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/eav.global.tsv
TE_consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
family_renames=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt

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
# Remove trailing tab and replace it with a newline
sed  's/\t$/\n/' column_names.txt> tmp && mv tmp column_names.txt
# Add on appropriately ordered column names to FINAL_RESULTS_LF_REF
cat column_names.txt FINAL_RESULTS_LF_ReF.txt >Full_Results.txt
# Add on appropriately ordered column names to FINAL_RESULTS.txt
cat te_names_ordered.txt FINAL_RESULTS.txt > tmp && mv tmp FINAL_RESULTS.txt
# Transpose the matrix, remove QX strains if still present
cat Full_Results.txt | awk '$1 != "QX2265" && $1 != "QX2266" {print $0}' > tmp && mv tmp Full_Results.txt
bash ${scripts_dir}/${transpose} Full_Results.txt T_Full_Results.txt
# Add information on whether the TE family is a dnatransposon, retrotransposon, or unknown
python ${scripts_dir}/${assign_cut_copy} $repbase_fasta $consensus_renamed ${data_dir}/all_nonredundant.txt
# remove lines with all zeros (no instance of that transposition event in any of the samples)
mv T_Full_Results.txt original_T_Full_Results.txt
cat original_T_Full_Results.txt | awk '{for (i=2;i<=NF;i++) {if($i !=0){print $0;break}}}' > T_Full_Results.txt
# calculate "NAs" for missing data
bash ${scripts_dir}/${NA} $samples
# generate kinship matrix for insertion calls
mkdir kinship
cd kinship
python ${scripts_dir}/${kin} $samples ${results_dir}/CtCp_all_nonredundant.txt
mv TE_matrix/kin_matrix_ins.txt .
#generate kinship matrix for reference/absence calls
python ${scripts_dir}/${kin_AF} $samples
#remove monomorphic reference calls (all strains have that transposon(or an NA))
cp kin_matrix_AF.txt with_monomorphic_calls_kin_matrix_AF.txt
cat kin_matrix_AF.txt | awk '$0~/\t0\t/||$0~/trait/ {print $0}' >tmp && mv tmp kin_matrix_AF.txt
# make sure column names between the matrices match
result=`diff <(head -n 1 kin_matrix_ins.txt) <(head -n 1 kin_matrix_AF.txt)`
if [$result eq '']
	then echo "matched"
	else echo "column names do not match, exiting..." && exit 1
fi
# merge kinship matrices:
cat kin_matrix_ins.txt > tmp && cat kin_matrix_AF.txt |sed 1d >>tmp && mv tmp kin_matrix_full.txt
#transpose individual matrices:
bash ${scripts_dir}/transpose_matrix.sh kin_matrix_ins.txt  T_kin_matrix_ins.txt 
bash ${scripts_dir}/transpose_matrix.sh kin_matrix_AF.txt  T_kin_matrix_AF.txt 
#calculate percentage of NAs at each insertion position
python ${scripts_dir}/${NA_pos}
#reduce positions for mappings (remove sites with a high amount of missing data)
python ${scripts_dir}/${reduce_ins}
bash ${scripts_dir}/transpose_matrix.sh kin_matrix_ins_reduced.txt  T_kin_matrix_ins_reduced.txt 
# make sure column names between the matrices match
result=`diff <(head -n 1 kin_matrix_ins_reduced.txt) <(head -n 1 kin_matrix_AF.txt)`
if [$result eq '']
	then echo "matched"
	else echo "column names do not match, exiting..." && exit 1
fi
# merge kinship matrices:
cat kin_matrix_ins_reduced.txt > tmp && cat kin_matrix_AF.txt |sed 1d >>tmp && mv tmp kin_matrix_full_reduced.txt



#calculate fraction totals
python ${scripts_dir}/${frac_matrix_ins} T_kin_matrix_ins.txt
python ${scripts_dir}/${frac_matrix_AF} T_kin_matrix_AF.txt
python ${scripts_dir}/${frac_matrix_ins_reduced} T_kin_matrix_ins_reduced.txt
#transpose matrices
bash ${scripts_dir}/transpose_matrix.sh kin_frac_matrix_ins.txt T_kin_frac_matrix_ins.txt
bash ${scripts_dir}/transpose_matrix.sh kin_frac_matrix_AF.txt T_kin_frac_matrix_AF.txt
bash ${scripts_dir}/transpose_matrix.sh kin_C_matrix_ins.txt T_kin_C_matrix_ins.txt
bash ${scripts_dir}/transpose_matrix.sh kin_C_matrix_AF.txt T_kin_C_matrix_AF.txt
bash ${scripts_dir}/transpose_matrix.sh kin_C_matrix_NAs.txt T_kin_C_matrix_NAs.txt
#transpose reduced matrices
bash ${scripts_dir}/transpose_matrix.sh kin_frac_matrix_ins_reduced.txt T_kin_frac_matrix_ins_reduced.txt
bash ${scripts_dir}/transpose_matrix.sh kin_C_matrix_ins_reduced.txt T_kin_C_matrix_ins_reduced.txt
bash ${scripts_dir}/transpose_matrix.sh kin_C_matrix_NAs_reduced.txt T_kin_C_matrix_NAs_reduced.txt

# make sure column names between the matrices match
result=`diff <(head -n 1 T_kin_frac_matrix_ins.txt) <(head -n 1 T_kin_frac_matrix_AF.txt)`
if [$result eq '']
	then echo "matched"
	else echo "column names do not match, exiting..." && exit 1
fi

# make sure column names between the matrices match
result=`diff <(head -n 1 T_kin_frac_matrix_ins_reduced.txt) <(head -n 1 T_kin_frac_matrix_AF.txt)` #AF matrix doesn't need to be reduced
if [$result eq '']
	then echo "matched"
	else echo "column names do not match, exiting..." && exit 1
fi
#merge matrices
cat T_kin_frac_matrix_ins.txt > tmp && cat T_kin_frac_matrix_AF.txt |sed 1d >>tmp && mv tmp T_kin_frac_matrix_full.txt
cat T_kin_C_matrix_ins.txt > tmp && cat T_kin_C_matrix_AF.txt |sed 1d >>tmp && mv tmp T_kin_C_matrix_full.txt
#merge reduced insertion matrix with original AF matrix
cat T_kin_frac_matrix_ins_reduced.txt > tmp && cat T_kin_frac_matrix_AF.txt |sed 1d >>tmp && mv tmp T_kin_frac_matrix_full_reduced.txt
cat T_kin_C_matrix_ins_reduced.txt > tmp && cat T_kin_C_matrix_AF.txt |sed 1d >>tmp && mv tmp T_kin_C_matrix_full_reduced.txt
# remove lines with all NA (no instance of that transposition event in any of the samples)
mv T_kin_frac_matrix_full.txt original_T_kin_frac_matrix_full.txt
cat original_T_kin_frac_matrix_full.txt | awk '{for (i=2;i<=NF;i++) {if($i !="NA"){print $0;break}}}' > T_kin_frac_matrix_full.txt
mv T_kin_C_matrix_full.txt original_T_kin_C_matrix_full.txt
cat original_T_kin_C_matrix_full.txt | awk '{for (i=2;i<=NF;i++) {if($i !="NA"){print $0;break}}}' > T_kin_C_matrix_full.txt
# remove lines with all NA  for reduced files
mv T_kin_frac_matrix_full_reduced.txt original_T_kin_frac_matrix_full_reduced.txt
cat original_T_kin_frac_matrix_full_reduced.txt | awk '{for (i=2;i<=NF;i++) {if($i !="NA"){print $0;break}}}' > T_kin_frac_matrix_full_reduced.txt
mv T_kin_C_matrix_full_reduced.txt original_T_kin_C_matrix_full_reduced.txt
cat original_T_kin_C_matrix_full_reduced.txt | awk '{for (i=2;i<=NF;i++) {if($i !="NA"){print $0;break}}}' > T_kin_C_matrix_full_reduced.txt
#transpose and calculate activity measurements
bash ${scripts_dir}/transpose_matrix.sh T_kin_C_matrix_full.txt kin_C_matrix_full.txt
bash ${scripts_dir}/transpose_matrix.sh T_kin_C_matrix_full_reduced.txt kin_C_matrix_full_reduced.txt
python ${scripts_dir}/${activity} kin_C_matrix_full.txt #only for non-redcued set
bash ${scripts_dir}/transpose_matrix.sh Full_Results_Activity.txt T_Full_Results_Activity.txt
# add in coverage trait
python ${scripts_dir}/${coverage} ../Full_Results.txt $bam_stats
python ${scripts_dir}/${fix_cov}

result=`diff <(head -n 1 T_kin_C_matrix_full.txt) <(head -n 1 T_coverage.txt )`
if [$result eq '']
	then echo "matched"
	else echo "column names do not match, exiting..." && exit 1
fi

result=`diff <(head -n 1 T_kin_C_matrix_full_reduced.txt) <(head -n 1 T_coverage.txt )`
if [$result eq '']
	then echo "matched"
	else echo "column names do not match, exiting..." && exit 1
fi
cat T_coverage.txt |sed '1d' >tmp
cat T_kin_C_matrix_full.txt tmp > all && mv all T_kin_C_matrix_full.txt
cat T_kin_C_matrix_full_reduced.txt tmp > all && mv all T_kin_C_matrix_full_reduced.txt
# copy files to final folder
cd ${results_dir}
mkdir final_results
cp ${results_dir}/kinship/T_Full_Results_Activity.txt ${results_dir}/kinship/T_kin_frac_matrix_full.txt ${results_dir}/kinship/T_kin_C_matrix_full.txt ${results_dir}/kinship/kin_matrix_full.txt ${results_dir}/kinship/kin_matrix_full_reduced.txt ${results_dir}/kinship/T_kin_C_matrix_full_reduced.txt final_results
mv CtCp_all_nonredundant.txt CtCp_redundant.txt # rename this file to avoid confustion
#generate CtCp file with assigned TE classes for both unreduced and reduced datasets
cd final_results
python ${scripts_dir}/generate_CtCp.py kin_matrix_full.txt
python ${scripts_dir}/generate_CtCp_reduced.py kin_matrix_full_reduced.txt
python ${scripts_dir}/${assign_cut_copy_R2} $repbase_fasta $consensus_renamed all_nonredundant.txt
python ${scripts_dir}/${assign_cut_copy_R2} $repbase_fasta $consensus_renamed all_nonredundant_reduced.txt
rm all_nonredundant.txt
rm all_nonredundant_reduced.txt

cat T_Full_Results_Activity.txt | awk '{for (i=2;i<=NF;i++) {if($i !="NA"){print $0;break}}}' > tmp && mv tmp T_Full_Results_Activity.txt 
cat kin_matrix_full.txt | awk '{for (i=2;i<=NF;i++) {if($i !="NA"){print $0;break}}}' > tmp && mv tmp kin_matrix_full.txt
cat kin_matrix_full_reduced.txt | awk '{for (i=2;i<=NF;i++) {if($i !="NA"){print $0;break}}}' > tmp && mv tmp kin_matrix_full_reduced.txt
#clip CtCp files to remove redundancies
cat CtCp_all_nonredundant.txt |sort -k1,1 -k2,2n -k3,3 -k4,4 -k6,6|awk '!x[$1,$2,$3,$4,$6]++' > CtCp_clipped.txt #LEFT OFF HERE
cat CtCp_clipped.txt| sort -k1,1 -k2,2n > tmp && mv tmp CtCp_clipped.txt
cat CtCp_clipped.txt |awk '{print $1"\tTE\t"$4"\t"$2+1"\t"$3+1"\t"$6"\t"$5"\tNA\t"$8}' > CtCp_clipped.gff 
#clip CtCp files to remove redundancies
cat CtCp_all_nonredundant_reduced.txt |sort -k1,1 -k2,2n -k3,3 -k4,4 -k6,6|awk '!x[$1,$2,$3,$4,$6]++' > CtCp_clipped_reduced.txt #LEFT OFF HERE
cat CtCp_clipped_reduced.txt| sort -k1,1 -k2,2n > tmp && mv tmp CtCp_clipped_reduced.txt
cat CtCp_clipped_reduced.txt |awk '{print $1"\tTE\t"$4"\t"$2+1"\t"$3+1"\t"$6"\t"$5"\tNA\t"$8}' > CtCp_clipped_reduced.gff 
#determine genomic features that TEs are located in 
bash ${scripts_dir}/${gene_interrupt}
#CER1 checks
python ${scripts_dir}/${cerfinder} #reduction doesn't affect ref/ab calls...doesn't matter which set used here
#find outliers for total ref,abs, ins, dna, retro, unknown
python ${scripts_dir}/${find_outliers}
#put input files for figure generation in new directory
python ${scripts_dir}/${count_classes} 
result=`diff <(head -n 1 T_kin_C_matrix_full.txt) <(head -n 1 total_classes.txt )`
if [$result eq '']
	then echo "matched"
	else echo "column names do not match, exiting..." && exit 1
fi
cat T_kin_C_matrix_full.txt > tmp && cat total_classes.txt |sed 1d >>tmp && mv tmp T_kin_C_matrix_full.txt


python ${scripts_dir}/${count_classes_reduced} #REDUCED
result=`diff <(head -n 1 T_kin_C_matrix_full_reduced.txt) <(head -n 1 total_classes_reduced.txt )`
if [$result eq '']
	then echo "matched"
	else echo "column names do not match, exiting..." && exit 1
fi
cat T_kin_C_matrix_full_reduced.txt > tmp && cat total_classes_reduced.txt |sed 1d >>tmp && mv tmp T_kin_C_matrix_full_reduced.txt

mkdir data_for_figures
#genes of interest
cd gene_interrupt
mkdir models
python ${scripts_dir}/${interest} vit-1 cdr-2  #gon-2 prg-1 unc-130 ceh-34 egl-34 pkd-2 --don't need this step anymore
cat essentiality_nonredundant_cds.txt |awk '$9 !="NA" {print $0}' > essentiality_nonredundant_cds_lethal.txt


# per strain NA pruning
cd ${results_dir}/final_results
cp T_kin_C_matrix_full_reduced.txt T_kin_C_matrix_full_reduced_unpruned.txt
python ${scripts_dir}/NA_prune.py 
cp pruned_data.txt T_kin_C_matrix_full_reduced.txt

# check which TE families predominate the total insertion count
python ${scripts_dir}/count_predominate_insTE.py > predominate_insertions.txt
python ${scripts_dir}/check_nums.py $repbase_fasta $consensus_renamed CtCp_all_nonredundant.txt


#simplify names (WBTransposon->WBT, remove trailing _CE)
#cat kin_matrix_full.txt | awk  '{if($1~/non-reference/){gsub("non-reference.*","NR",$1);gsub("_CE_NR$","_NR",$1)}else{$1=$1"_R";gsub("_CE_R$","_R",$1)}; print $0;}' |sed 's/WBTransposon/WBT/' > tmp && mv tmp kin_matrix_full.txt 
cat kin_matrix_full.txt  | awk  '{if($1~/non-reference/){gsub("non-reference.*","NR",$1);gsub("_CE_NR$","_NR",$1)}else{if($1!~/trait/){$1=$1"_R";gsub("_CE_R$","_R",$1)}}; print $0;}' |sed 's/WBTransposon/WBT/' > tmp && mv tmp kin_matrix_full.txt 
cat T_kin_C_matrix_full.txt |sed 's/WBTransposon/WBT/'|sed 's/_CE_C/_C/' >tmp && mv tmp T_kin_C_matrix_full.txt
cat T_kin_C_matrix_full_reduced.txt |sed 's/WBTransposon/WBT/'|sed 's/_CE_C/_C/'|awk -v OFS="\t" '{if($1~/_CE$/) {gsub("_CE$","",$1)};print $0}' > tmp && mv tmp T_kin_C_matrix_full_reduced.txt

#cat T_kin_C_matrix_full_reduced.txt |sed 's/WBTransposon/WBT/'|sed 's/_CE_C/_C/'|awk -v OFS="\t" '{if($1~/_CE$/) {gsub("_CE$","",$1)};print $0}'>test.txt



cat CtCp_all_nonredundant.txt |sed 's/WBTransposon/WBT/'|sed 's/_CE_/_/' > tmp && mv tmp CtCp_all_nonredundant.txt
cat CtCp_all_nonredundant_reduced.txt |sed 's/WBTransposon/WBT/'|sed 's/_CE_/_/' > tmp && mv tmp CtCp_all_nonredundant_reduced.txt


#assign id numbers to traits
python ${scripts_dir}/${id}
python ${scripts_dir}/${id_reduced}

#get supplementary fasta table
python ${scripts_dir}/${fasta_supp} $TE_consensus $family_renames

cd ${results_dir}/final_results/data_for_figures
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/coverage_and_te_counts.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/T_kin_C_matrix_ins_reduced.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/NA_counts_at_positions.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/T_kin_C_matrix_NAs_reduced.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/essentiality_nonredundant_GO.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/essentiality_nonredundant_cds.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/kin_matrix_full.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/contradictory_calls.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full_reduced.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/CtCp_all_nonredundant.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/CtCp_all_nonredundant_reduced.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/cer_comparison.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/key_T_kin_C_matrix_full_id.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full_id_reduced.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/key_T_kin_C_matrix_full_id_reduced.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/key_kin_matrix_full_id.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/key_T_Full_Results_Activity_id.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/outliers_fam_tot.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/outliers_families_pruned.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/TE_seqs.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/homologs.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/paragraphs.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/grantham_scores.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/FINAL_RESULTS_LF.txt_IGNORE .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/OD/OD_FINAL_RESULTS_LF.txt .





echo "DONE"
