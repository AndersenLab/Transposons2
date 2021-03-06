#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --nodelist=node4
#SBATCH --mem=2000

python_version=/exports/people/andersenlab/kml436/python2/bin/python
bam_surgeon=/lscr2/andersenlab/kml436/git_repos2/bamsurgeon
reference=/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa
run_ID=run_${1}
bam=/lscr2/andersenlab/dec211/RUN/v2_snpset/bam/N2.bam
TE_consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
processors=18
TTR=/lscr2/andersenlab/kml436/git_repos2/mcclintock/
#HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_tes.gff ###
HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names_alias.bed
original_ref_pos=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names.gff #actually bed file
element_list=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round19/fasta_list.txt ###
location_list=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round19/bed_list.txt ###
minimal_Distance_to_count=1000
minimal_supporting_reads=3
minimal_supporting_individuals=1
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt
fake_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths_plusTC8.txt ####
high_coverage_regions=/lscr2/andersenlab/kml436/sv_files/N2_high_coverage_regions
existing_transposons=/lscr2/andersenlab/kml436/ID_transposable_element.bed
family_renames=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt 
all_WB_seqs=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/Wb-TC8.fasta
ref_tes=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_family_names.gff
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 
bam_name=$(basename "$bam" .bam)
echo "bam name is ${bam_name}"
mkdir ${run_ID}_${bam_name}
dir1=`pwd`
cd ${run_ID}_${bam_name}
dir=`pwd`
##############################################################
##create error log file for everything?
##randomsites--add options here?

find . -name "*unpair*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*insertion*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*align*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*uniq*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*clipped*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*filter_results*" -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*final_results*" -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
mkdir ${run_ID}_${bam_name}_non_filter_results/
mkdir ${run_ID}_${bam_name}_filter_results/

#mv ${run_ID}_${bam_name}.sorted.bam copy_${run_ID}_${bam_name}.sorted.bam
# ${run_ID}_${bam_name}.sorted.bam.bai copy_${run_ID}_${bam_name}.sorted.bam.bai
#quarter of the reads
#java -jar /opt/picard-tools/DownsampleSam.jar P=.25 I=copy_${run_ID}_${bam_name}.sorted.bam O=testing.bam
#/opt/samtools/x64/samtools/samtools sort -o -@ 8 testing.bam out > ${run_ID}_${bam_name}.sorted.bam
#/opt/samtools/x64/samtools/samtools index ${run_ID}_${bam_name}.sorted.bam


#java -jar /opt/picard-tools/DownsampleSam.jar P=.5 I=copy_${run_ID}_${bam_name}.sorted.bam O=testing.bam
#/opt/samtools/x64/samtools/samtools sort -o -@ 8 testing.bam out > ${run_ID}_${bam_name}.sorted.bam
#/opt/samtools/x64/samtools/samtools index ${run_ID}_${bam_name}.sorted.bam

##############################################################
echo "Running TEMP..."
bash ${TTR}/TEMP/scripts/TEMP_Insertion.sh -i ${run_ID}_${bam_name}.sorted.bam -s ${TTR}/TEMP/scripts/ -x 30 -r $TE_consensus -t $original_ref_pos -m 1 -c 20 &> ${run_ID}_TEMP_log.txt
echo "TEMP finished...altering output..."
sed '1d' "${run_ID}_${bam_name}.insertion.refined.bp.summary" | awk '{if ($4 == "sense" || $4 == "antisense"); else print $0}' | awk '{ printf $1"\t"$9"\t"$11"\t"$7"\t"$4"_non-reference_"sample"_temp_\t0\t"; if ( $5 == "sense" ) printf "+"; else printf "-"; print "\t"$6"\t"$10"\t"$12"\t"$8}' > "${run_ID}_${bam_name}_temp_presort_raw.txt"
bedtools sort -i "${run_ID}_${bam_name}_temp_presort_raw.txt" > "${run_ID}_${bam_name}_temp_sorted_raw.txt"
awk '{ printf $1"\t"$2"\t"$3"\t"$4"\t"$5; if ($9 > 0 && $10 > 0) printf "sr_"; else printf "rp_"; print NR"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' "${run_ID}_${bam_name}_temp_sorted_raw.txt" > "${run_ID}_${bam_name}_temp_sorted_precut_raw.txt"
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7"\t"$11}' "${run_ID}_${bam_name}_temp_sorted_precut_raw.txt" >> "${run_ID}_${bam_name}_temp_raw.bed"
mv ${run_ID}_${bam_name}_temp_raw.bed ${run_ID}_${bam_name}_non_filter_results/
echo "Filtering TEMP output......"
# Filter for insertions with support for both ends
awk '{if ($8 == "1p1") print $0}' "${run_ID}_${bam_name}_temp_sorted_precut_raw.txt" > "${run_ID}_${bam_name}_temp_sorted_redundant.txt"
cut -f1-3,5-7 "${run_ID}_${bam_name}_temp_sorted_redundant.txt" >> "${run_ID}_${bam_name}_temp_redundant.bed"
# Filter out redundant insertions
sort -k1,3 -k4rn "${run_ID}_${bam_name}_temp_sorted_redundant.txt" | sort -u -k1,3 | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7"\t"$11}' > "${run_ID}_${bam_name}_tmp"
bedtools sort -i "${run_ID}_${bam_name}_tmp">> "${run_ID}_${bam_name}_temp_nonredundant.bed"
mv ${run_ID}_${bam_name}_temp_nonredundant.bed ${run_ID}_${bam_name}_filter_results/
cd ${run_ID}_${bam_name}_non_filter_results
mv ${run_ID}_${bam_name}_temp_raw.bed  ${run_ID}_${bam_name}_temp_nonredundant.bed
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_temp_nonredundant.bed TMP
mv ${run_ID}_${bam_name}_temp_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy.bed
mv TMP ${run_ID}_${bam_name}_temp_nonredundant.bed
cat ${run_ID}_${bam_name}_temp_nonredundant.bed| sed s'/chrX/X/g'| sed '/^track/d' | awk '$4 ~/non-reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${bam_name}_temp && mv ${bam_name}_temp ${run_ID}_${bam_name}_temp_nonredundant.bed  

cd ..
cd ${run_ID}_${bam_name}_filter_results
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_temp_nonredundant.bed TMP
mv ${run_ID}_${bam_name}_temp_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy.bed
mv TMP ${run_ID}_${bam_name}_temp_nonredundant.bed
cat ${run_ID}_${bam_name}_temp_nonredundant.bed| sed s'/chrX/X/g'| sed '/^track/d' | awk '$4 ~/non-reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${bam_name}_temp && mv ${bam_name}_temp ${run_ID}_${bam_name}_temp_nonredundant.bed  

cd ..
##############################################################


##############################################################

echo "Comparing outputs to simulated structural variants......"
#python /exports/people/andersenlab/kml436/scripts/create_bam_surgeon_log.py ${run_ID}_${bam_name}_bamsurgeon.log
#bedtools sort -i ${run_ID}_${bam_name}_final_positions.bed> tmpSUM && mv tmpSUM ${run_ID}_${bam_name}_final_positions.bed
#cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_final_positions_WB_copy.bed 
#python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2.py $family_renames ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_final_positions_TMP.bed
#mv ${run_ID}_${bam_name}_final_positions_TMP.bed ${run_ID}_${bam_name}_final_positions.bed
cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_non_filter_results/
cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_filter_results/

mkdir final_results

cd ${run_ID}_${bam_name}_non_filter_results/
mkdir BEDCOMPARE_files
${python_version} ${bam_surgeon}/bed_compare9_redo_EDIT_round20.py 20 ${dir}/${run_ID}_${bam_name}_non_filter_results/ $fake_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} $consensus_renamed >& ${run_ID}_${bam_name}_bedcompare.log
cd intermediates/
cat FAMILY_TFPN* > FAMILY_TFPNs_NF
mv FAMILY_TFPNs_NF ../../final_results
cd ..
mv BEDCOMPARE_${run_ID}_${bam_name}.txt BEDCOMPARE_FAMILY_${run_ID}_${bam_name}.txt BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt BEDCOMPARE_files
cd BEDCOMPARE_files
cat BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt | sed 's/redundant\.bed/redundant\.bed_NF/g' > tmp_non_filter_file.bed
cd ../..

cd ${run_ID}_${bam_name}_filter_results/
mkdir BEDCOMPARE_files
${python_version} ${bam_surgeon}/bed_compare9_redo_EDIT_round20.py 20 ${dir}/${run_ID}_${bam_name}_filter_results/ $fake_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} $consensus_renamed >& ${run_ID}_${bam_name}_bedcompare.log
cd intermediates/
cat FAMILY_TFPN* > FAMILY_TFPNs_F
mv FAMILY_TFPNs_F ../../final_results
cd ..
mv BEDCOMPARE_${run_ID}_${bam_name}.txt BEDCOMPARE_FAMILY_${run_ID}_${bam_name}.txt BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt BEDCOMPARE_files
cd BEDCOMPARE_files
cat BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt | sed 1d | sed 's/redundant\.bed/redundant\.bed_F/g' > tmp_filter_file.bed
cd ../..

cat ${run_ID}_${bam_name}_non_filter_results/BEDCOMPARE_files/tmp_non_filter_file.bed ${run_ID}_${bam_name}_filter_results/BEDCOMPARE_files/tmp_filter_file.bed > BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt
mv BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt final_results

echo "DONE"