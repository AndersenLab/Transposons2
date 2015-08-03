#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=18000

# this script runs round18 of the transposon simulations/detection
# USE: WB_SIM_round18_TransposonsBS_N2_SB_TTR.sh

python_version=/exports/people/andersenlab/kml436/python2/bin/python
bam_surgeon=/lscr2/andersenlab/kml436/git_repos2/bamsurgeon
reference=/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa
run_ID=run_${1}
bam=/lscr2/andersenlab/dec211/RUN/v2_snpset/bam/N2.bam
processors=18
TTR=/lscr2/andersenlab/kml436/git_repos2/mcclintock/
#HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round2_WB_tes.gff
HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/WB_pos_element_names_alias.bed
minimal_Distance_to_count=1000
minimal_supporting_reads=3
minimal_supporting_individuals=1
high_coverage_regions=/lscr2/andersenlab/kml436/sv_files/N2_high_coverage_regions
existing_transposons=/lscr2/andersenlab/kml436/ID_transposable_element.bed

all_WB_seqs=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/Wb-TC8.fasta
ref_tes=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/WB_pos_family_names.gff



TE_consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/round2_consensus_set2.fasta
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/AB-PR/consensus_wTC8.fasta 
#consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/LENGTHS/consensus.fasta
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/AB-PR/fake_lengths.txt
#TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/LENGTHS/lengths.txt
WB_elements=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/Wb-TC8.fasta #these are the plus strands without TC8
family_renames=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round2_WB_familes_set2.txt
original_ref_pos=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/WB_pos_element_names.gff
bam_name=$(basename "$bam" .bam)
echo "bam name is ${bam_name}"
mkdir ${run_ID}_${bam_name}
dir=`pwd`
cd ${run_ID}_${bam_name}


##############################################################
##create error log file for everything?
##randomsites--add options here?


##spike ins (using altered csripts)y

find . -name "*unpair*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*insertion*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*align*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*uniq*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*clipped*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*filter_results*" -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*final_results*" -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*temp*" -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
mkdir ${run_ID}_${bam_name}_non_filter_results/
mkdir ${run_ID}_${bam_name}_filter_results/
##############################################################
#############################################################
echo "Creating sam file..."
mkdir ${run_ID}_${bam_name}_sam_file
samtools view ${run_ID}_${bam_name}.sorted.bam |sort --temporary-directory=${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_sam_file > ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_sam_file/${run_ID}_${bam_name}.sorted.sam
##in  /lscr2/andersenlab/kml436/git_repos2/mcclintock/TE-locate/
echo "running TELOCATE..."

mkdir TELOCATE
cd ${TTR}/TE-locate/
perl TE_locate.pl 2 ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_sam_file/ $HL_gff $reference ${dir}/${run_ID}_${bam_name}/TELOCATE/TEL $minimal_Distance_to_count $minimal_supporting_reads $minimal_supporting_individuals
cd ${dir}/${run_ID}_${bam_name}
echo "TELOCATE finished...altering output..."
##LEFT OFF HERE
sed -e '1,2d' "${dir}/${run_ID}_${bam_name}/TELOCATE/TEL_${minimal_Distance_to_count}_reads${minimal_supporting_reads}_acc${minimal_supporting_individuals}.info" > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile"
awk -F'[\t/]' '{printf $1"\t"; if($17=="old") printf $2-1"\t"$2"\t"$8"\t"$5"_reference_telocate_rp_\t0"; else if($17=="new") printf $2-1"\t"$2"\t"$8"\t"$5"_non-reference_telocate_rp_\t0"; if($14=="parallel") print "\t+"; else if($14 =="uncertain") print "\t."; else print "\t-";}' "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile" > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_presort.txt"
bedtools sort -i "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_presort.txt" > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_sorted_redundant.txt"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7}' "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_sorted_redundant.txt" > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_counted_redundant.txt"
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_counted_redundant.txt" >> "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed"
cat ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed | awk '$4~/_reference/ {print $0}' > tmp && mv tmp ${run_ID}_tmpfile_telocate_redundant.bed 
mv ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed ${run_ID}_${bam_name}_non_filter_results/
echo "STEP1"
cd ${run_ID}_${bam_name}_non_filter_results/
#rename TELOCATE NON FILTER results (from WB element names to consensus names)....won't actualy end up using the "unfiltered" results...same as "filtered"
mv ${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed ${run_ID}_${bam_name}_telocate_nonredundant.bed 
python /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_telocate_nonredundant.bed  TMP
mv ${run_ID}_${bam_name}_telocate_nonredundant.bed  ${run_ID}_nonredundant_org_copy.bed
mv TMP ${run_ID}_${bam_name}_telocate_nonredundant.bed 
cat ${run_ID}_${bam_name}_telocate_nonredundant.bed| awk '$4 ~/_reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${run_ID}_temp && mv ${run_ID}_temp ${run_ID}_${bam_name}_telocate_nonredundant.bed 

cd ..
echo "STEP2"
echo "Filtering TELOCATE output......"
sort -k1,3 -k4rn "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_counted_redundant.txt"| sort -u -k1,3 | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfileNR"
echo "STEP3"
bedtools sort -i "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfileNR" >> "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed"
echo "STEP4"
cat ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed| awk '$4~/_reference/ {print $0}' > tmp && mv tmp ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed
mv ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed ${run_ID}_${bam_name}_filter_results/
######
echo "STEP5"
echo "STEP6"
echo "STEP7"
#rename TELOCATE results (from WB element names to consensus names)
cd ${run_ID}_${bam_name}_filter_results/
python /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_telocate_nonredundant.bed  TMP
mv ${run_ID}_${bam_name}_telocate_nonredundant.bed  ${run_ID}_nonredundant_org_copy.bed
mv TMP ${run_ID}_${bam_name}_telocate_nonredundant.bed

#find closest matches to reference TEs within 1000 base pairs
closestBed -a $ref_tes -b ${run_ID}_${bam_name}_telocate_nonredundant.bed -d -t all | awk '$13<=1000 {print $0}' > ${run_ID}_closest.txt
python /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/remove_telocate_redundancies.py ${run_ID}_closest.txt
python /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/process_telocate_output.py closest_prepped.txt
cat processed_calls.txt| sort -k1,1 -k2,2n > ${run_ID}_${bam_name}_telocate_nonredundant.bed
cat ${run_ID}_${bam_name}_telocate_nonredundant.bed| awk '$4 ~/_reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}' > ${run_ID}_temp && mv ${run_ID}_temp ${run_ID}_${bam_name}_telocate_nonredundant.bed
cd ..

echo "STEP8"
# get the  positions of the simulated reference TEs
cat $ref_tes | awk '$5=="R" && $4 !="TC8" {print $0}' > tmp && mv tmp ref_pos.bed
cp ref_pos.bed ${run_ID}_${bam_name}_non_filter_results
cp ref_pos.bed ${run_ID}_${bam_name}_filter_results/
echo "STEP9"
#FEED JUST REFERENCE to TELOCATE
echo "Comparing outputs to simulated structural variants......"
mkdir final_results_RECALCULATED
cd ${run_ID}_${bam_name}_filter_results/
mkdir BEDCOMPARE_files_RECALCULATED
python /lscr2/andersenlab/kml436/git_repos2/bamsurgeon/bed_compare9_redo_EDIT_RECALCULATE_absence.py 20 ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_filter_results/ $TE_lengths ref_pos.bed ${run_ID}_${bam_name} $consensus_renamed >& ${run_ID}_${bam_name}_bedcompare.log
cd intermediates_RECALCULATED/
cat FAMILY_TFPN* > FAMILY_TFPNs_telocate
mv FAMILY_TFPNs_telocate ../../final_results_RECALCULATED
cd ..
mv BEDCOMPARE_${run_ID}_${bam_name}.txt BEDCOMPARE_FAMILY_${run_ID}_${bam_name}.txt BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt BEDCOMPARE_files_RECALCULATED
cd BEDCOMPARE_files_RECALCULATED
cat BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt |sed 1d | sed 's/redundant\.bed/redundant\.bed_F/g' > tmp_filter_file.bed
cd ../..
echo "STEP10"

#cd ${run_ID}_${bam_name}_non_filter_results/
#mkdir BEDCOMPARE_files_RECALCULATED
#python /lscr2/andersenlab/kml436/git_repos2/bamsurgeon/bed_compare9_redo_EDIT_RECALCULATE_absence.py 20 ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_non_filter_results/ $TE_lengths ref_pos.bed ${run_ID}_${bam_name} $consensus_renamed >& ${run_ID}_${bam_name}_bedcompare.log
#cd intermediates_RECALCULATED/
#cat FAMILY_TFPN* > FAMILY_TFPNs_NF
#mv FAMILY_TFPNs_NF ../../final_results_RECALCULATED
#cd ..
#mv BEDCOMPARE_${run_ID}_${bam_name}.txt BEDCOMPARE_FAMILY_${run_ID}_${bam_name}.txt BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt BEDCOMPARE_files_RECALCULATED
#cd BEDCOMPARE_files_RECALCULATED
#cat BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt |sed 1d | sed 's/redundant\.bed/redundant\.bed_NF/g' > tmp_non_filter_file.bed
#cd ../..

#########################

echo "STEP11"
cat ${run_ID}_${bam_name}_filter_results/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed > BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt
#${run_ID}_${bam_name}_non_filter_results/BEDCOMPARE_files_RECALCULATED/tmp_non_filter_file.bed
mv BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt final_results_RECALCULATED
cd final_results_RECALCULATED
mv FAMILY_TFPNs_telocate FAMILY_TFPNs_F
echo "Done"

