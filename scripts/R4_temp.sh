#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=18000
#SBATCH --exclude=node[6,7,8]

# this script runs round4 of the transposon simulations/detection
# USE: WB_SIM_round4_TransposonsBS_N2_SB_TTR.sh

python_version=/exports/people/andersenlab/kml436/python2/bin/python
bam_surgeon=/lscr2/andersenlab/kml436/git_repos2/bamsurgeon
reference=/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa
run_ID=run_${1}
bam=/lscr2/andersenlab/dec211/RUN/v2_snpset/bam/N2.bam
TE_consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round4/round2_consensus_set2_round4.fasta
#TE_consensus_alternative=XXXXXXXXXXXX
##ALTERNATE: TE_consensus=/lscr2/andersenlab/kml436/common.fasta
processors=18
TTR=/lscr2/andersenlab/kml436/git_repos2/mcclintock/
element_list=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round4/fasta_list.txt
##ALTERNATE: element_list=/lscr2/andersenlab/kml436/sv_files/WB_REPB_elementlist.gff
location_list=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round4/bed_list.txt
##ALTERNATE: location_list=/lscr2/andersenlab/kml436/sv_files/WB_REPB_locationlist.gff
HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round4/WB_pos_element_names_alias_round4.bed
HL_gff_T=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round4/telcoate_WB_pos_element_names_alias_round4.bed
temp_ref=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round4/temp_ref_pos.bed
##ALTERNATE: HL_gff=/lscr2/andersenlab/kml436/sv_files/WB_REPB_ID_transposable_element_HL.gff
minimal_Distance_to_count=1000
minimal_supporting_reads=3
minimal_supporting_individuals=1
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt
high_coverage_regions=/lscr2/andersenlab/kml436/sv_files/N2_high_coverage_regions
existing_transposons=/lscr2/andersenlab/kml436/ID_transposable_element.bed
family_renames=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt
all_WB_seqs=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/Wb-TC8.fasta
bam_name=$(basename "$bam" .bam)
echo "bam name is ${bam_name}"
mkdir ${run_ID}_${bam_name}
cd ${run_ID}_${bam_name}
dir=`pwd`
##############################################################

echo "Running TEMP..."
bash ${TTR}/TEMP/scripts/TEMP_Insertion.sh -i ${run_ID}_${bam_name}.sorted.bam -s ${TTR}/TEMP/scripts/ -x 30 -r $TE_consensus -t $HL_gff -m 1 -c 20 &> ${run_ID}_TEMP_log.txt
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



mkdir final_results

cd ${run_ID}_${bam_name}_non_filter_results/
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_temp_nonredundant.bed TMP1
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_telocate_nonredundant.bed TMP2
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_retroseq_nonredundant.bed TMP3
mv ${run_ID}_${bam_name}_temp_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy1.bed
mv ${run_ID}_${bam_name}_telocate_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy2.bed
mv ${run_ID}_${bam_name}_retroseq_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy3.bed
mv TMP1 ${run_ID}_${bam_name}_temp_nonredundant.bed
mv TMP2 ${run_ID}_${bam_name}_telocate_nonredundant.bed
mv TMP3 ${run_ID}_${bam_name}_retroseq_nonredundant.bed
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
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_temp_nonredundant.bed TMP1
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_telocate_nonredundant.bed TMP2
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_retroseq_nonredundant.bed TMP3
mv ${run_ID}_${bam_name}_temp_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy1.bed
mv ${run_ID}_${bam_name}_telocate_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy2.bed
mv ${run_ID}_${bam_name}_retroseq_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy3.bed
mv TMP1 ${run_ID}_${bam_name}_temp_nonredundant.bed
mv TMP2 ${run_ID}_${bam_name}_telocate_nonredundant.bed
mv TMP3 ${run_ID}_${bam_name}_retroseq_nonredundant.bed
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