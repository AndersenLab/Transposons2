#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=18000

# this script runs round8 of the transposon simulations/detection
# USE: WB_SIM_round8_TransposonsBS_N2_SB_TTR.sh

python_version=/exports/people/andersenlab/kml436/python2/bin/python
bam_surgeon=/lscr2/andersenlab/kml436/git_repos2/bamsurgeon
reference=/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa
run_ID=run_${1}
bam=/lscr2/andersenlab/dec211/RUN/v2_snpset/bam/N2.bam
TE_consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/round2_consensus_fasta.fasta
#TE_consensus_alternative=XXXXXXXXXXXX
##ALTERNATE: TE_consensus=/lscr2/andersenlab/kml436/common.fasta
processors=18
TTR=/lscr2/andersenlab/kml436/git_repos2/mcclintock/
element_list=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/fasta_list.txt
##ALTERNATE: element_list=/lscr2/andersenlab/kml436/sv_files/WB_REPB_elementlist.gff
location_list=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/bed_list.txt
##ALTERNATE: location_list=/lscr2/andersenlab/kml436/sv_files/WB_REPB_locationlist.gff
HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round2_WB_tes.gff
##ALTERNATE: HL_gff=/lscr2/andersenlab/kml436/sv_files/WB_REPB_ID_transposable_element_HL.gff
minimal_Distance_to_count=1000
minimal_supporting_reads=3
minimal_supporting_individuals=1
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/lengths.txt
high_coverage_regions=/lscr2/andersenlab/kml436/sv_files/N2_high_coverage_regions
existing_transposons=/lscr2/andersenlab/kml436/ID_transposable_element.bed
family_renames=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round2_WB_familes.txt
all_WB_seqs=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/WB_all_seqs.fasta
bam_name=$(basename "$bam" .bam)
echo "bam name is ${bam_name}"
mkdir ${run_ID}_${bam_name}
cd ${run_ID}_${bam_name}
dir=`pwd`
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
mkdir ${run_ID}_${bam_name}_non_filter_results/
mkdir ${run_ID}_${bam_name}_filter_results/
##############################################################
echo "Running TEMP..."
bash ${TTR}/TEMP/scripts/TEMP_Insertion.sh -i ${run_ID}_${bam_name}.sorted.bam -s ${TTR}/TEMP/scripts/ -x 30 -r $TE_consensus -t $HL_gff -m 18 -c 20 &> ${run_ID}_TEMP_log.txt
echo "TEMP finished...altering output..."
sed '1d' "${run_ID}_${bam_name}.insertion.refined.bp.summary" | awk '{if ($4 == "sense" || $4 == "antisense"); else print $0}' | awk '{ printf $1"\t"$9"\t"$11"\t"$7"\t"$4"_non-reference_"sample"_temp_\t0\t"; if ( $5 == "sense" ) printf "+"; else printf "-"; print "\t"$6"\t"$10"\t"$12}' > "${run_ID}_${bam_name}_temp_presort_raw.txt"
bedtools sort -i "${run_ID}_${bam_name}_temp_presort_raw.txt" > "${run_ID}_${bam_name}_temp_sorted_raw.txt"
awk '{ printf $1"\t"$2"\t"$3"\t"$4"\t"$5; if ($9 > 0 && $10 > 0) printf "sr_"; else printf "rp_"; print NR"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' "${run_ID}_${bam_name}_temp_sorted_raw.txt" > "${run_ID}_${bam_name}_temp_sorted_precut_raw.txt"
cut -f1-3,5-7 "${run_ID}_${bam_name}_temp_sorted_precut_raw.txt" >> "${run_ID}_${bam_name}_temp_raw.bed"
mv ${run_ID}_${bam_name}_temp_raw.bed ${run_ID}_${bam_name}_non_filter_results/
echo "Filtering TEMP output......"
# Filter for insertions with support for both ends
awk '{if ($8 == "1p1") print $0}' "${run_ID}_${bam_name}_temp_sorted_precut_raw.txt" > "${run_ID}_${bam_name}_temp_sorted_redundant.txt"
cut -f1-3,5-7 "${run_ID}_${bam_name}_temp_sorted_redundant.txt" >> "${run_ID}_${bam_name}_temp_redundant.bed"
# Filter out redundant insertions
sort -k1,3 -k4rn "${run_ID}_${bam_name}_temp_sorted_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > "${run_ID}_${bam_name}_tmp"
bedtools sort -i "${run_ID}_${bam_name}_tmp">> "${run_ID}_${bam_name}_temp_nonredundant.bed"
mv ${run_ID}_${bam_name}_temp_nonredundant.bed ${run_ID}_${bam_name}_filter_results/

##############################################################
echo "running RETROSEQ-discovery..."
##Discovery
perl ${TTR}/RetroSeq/bin/retroseq.pl -discover -id 80 -bam ${run_ID}_${bam_name}.sorted.bam -eref $element_list -align -output ${run_ID}_${bam_name}_align.discovery
echo "running RETROSEQ-calling..."
##Calling
perl ${TTR}/RetroSeq/bin/retroseq.pl -call -bam ${run_ID}_${bam_name}.sorted.bam -input ${run_ID}_${bam_name}_align.discovery -filter $location_list -ref $reference -output ${run_ID}_${bam_name}_align.calling -orientate yes
echo "RETROSEQ finished...altering output..."
awk '$1!~/#/{print $0}' "${run_ID}_${bam_name}_align.calling.PE.vcf" > "Retro/tmp"
awk -F'[=,\t:]' '{print $1"\t"$11"\t"$12"\t"$6"\t"$10"_non-reference_retroseq_rp_\t0\t.\t"$21}' "Retro/tmp" > "${run_ID}_${bam_name}_align_retroseq_presort.txt"
bedtools sort -i "${run_ID}_${bam_name}_align_retroseq_presort.txt" > "${run_ID}_${bam_name}_align_retroseq_sorted_raw.txt"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7"\t"$8}' "${run_ID}_${bam_name}_align_retroseq_sorted_raw.txt" > "${run_ID}_${bam_name}_align_retroseq_sorted_counted_raw.txt"
cut -f1-3,5-7 "${run_ID}_${bam_name}_align_retroseq_sorted_counted_raw.txt" >> "${run_ID}_${bam_name}_align_retroseq_raw.bed"
mv ${run_ID}_${bam_name}_align_retroseq_raw.bed ${run_ID}_${bam_name}_non_filter_results/
## Filter results
echo "Filtering RETROSEQ output......"
awk '{if ($8 >= 6) print $0}' "${run_ID}_${bam_name}_align_retroseq_sorted_counted_raw.txt" > "${run_ID}_${bam_name}_align_retroseq_redundant.txt"
cut -f1-3,5-7 "${run_ID}_${bam_name}_align_retroseq_redundant.txt" >> "${run_ID}_${bam_name}_align_retroseq_redundant.bed"
## Filter for redundant predictions
sort -k1,3 -k4rn "${run_ID}_${bam_name}_align_retroseq_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > Retro/temp
bedtools sort -i Retro/temp >> "${run_ID}_${bam_name}_align_retroseq_nonredundant.bed"
mv ${run_ID}_${bam_name}_align_retroseq_nonredundant.bed ${run_ID}_${bam_name}_filter_results/

cd ${run_ID}_${bam_name}_non_filter_results
mv ${run_ID}_${bam_name}_temp_raw.bed  ${run_ID}_${bam_name}_temp_nonredundant.bed
mv ${run_ID}_${bam_name}_align_retroseq_raw.bed ${run_ID}_${bam_name}_retroseq_nonredundant.bed
cd ..
cd ${run_ID}_${bam_name}_filter_results
mv ${run_ID}_${bam_name}_align_retroseq_nonredundant.bed ${run_ID}_${bam_name}_retroseq_nonredundant.bed
cd ..

#########################
echo "Comparing outputs to simulated structural variants......"
#python /exports/people/andersenlab/kml436/scripts/create_bam_surgeon_log.py ${run_ID}_${bam_name}_bamsurgeon.log
#bedtools sort -i ${run_ID}_${bam_name}_final_positions.bed> tmpSUM && mv tmpSUM ${run_ID}_${bam_name}_final_positions.bed
#cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_final_positions_WB_copy.bed 
#python /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/correct_names_2.py $family_renames ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_final_positions_TMP.bed
#mv ${run_ID}_${bam_name}_final_positions_TMP.bed ${run_ID}_${bam_name}_final_positions.bed
cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_non_filter_results/
cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_filter_results/

mkdir final_results

cd ${run_ID}_${bam_name}_non_filter_results/
mkdir BEDCOMPARE_files
${python_version} ${bam_surgeon}/bed_compare9_redo_EDIT.py 20 ${dir}/${run_ID}_${bam_name}_non_filter_results/ $TE_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} $TE_consensus >& ${run_ID}_${bam_name}_bedcompare.log
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
${python_version} ${bam_surgeon}/bed_compare9_redo_EDIT.py 20 ${dir}/${run_ID}_${bam_name}_filter_results/ $TE_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} $TE_consensus >& ${run_ID}_${bam_name}_bedcompare.log
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

