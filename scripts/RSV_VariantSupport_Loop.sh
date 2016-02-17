#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

# this script filters TE caller output for varying amounts of read support (1-50 reads) and outputs the resulting TPR and FDR information files for RSVSIM
# can change read support, so double check what it is set to
# USE: ROUND18_ReadSupport_Loop.sh <run_ID>

consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 
#consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/consensus.fasta
#TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons22/files/SET2/AB-PR/fake_lengths.txt
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt
WB_elements=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/Wb-TC8.fasta #these are the plus strands without TC8
family_renames=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt


A_N2_1=/lscr2/andersenlab/kml436/fasta_subset/Rockman096-R096.2-N2-paired-1.fq
A_N2_2=/lscr2/andersenlab/kml436/fasta_subset/Rockman096-R096.2-N2-paired-2.fq
B_N2_1=/lscr2/andersenlab/kml436/fasta_subset/Princeton-P-N2_CGC-paired-1.fq
B_N2_2=/lscr2/andersenlab/kml436/fasta_subset/Princeton-P-N2_CGC-paired-2.fq
C_N2_1=/lscr2/andersenlab/kml436/fasta_subset/Uchicago-L001-N2Baer-paired-1.fq
C_N2_2=/lscr2/andersenlab/kml436/fasta_subset/Uchicago-L001-N2Baer-paired-2.fq

TEMP_scripts=/lscr2/andersenlab/kml436/git_repos2/mcclintock/TEMP/scripts/
absence=/lscr2/andersenlab/kml436/git_repos2/mcclintock/TEMP/scripts/TEMP_Absence.sh
faToTwoBit_dir=/exports/people/andersenlab/kml436
original_ref_pos=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names.gff
#original_ref_pos=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_consensus_RSV_sim.bed 
TTR=/lscr2/andersenlab/kml436/git_repos2/mcclintock/
minimal_Distance_to_count=1000
minimal_supporting_reads=3
minimal_supporting_individuals=1
run_ID=${1}
mkdir run_${run_ID}
cd run_${run_ID}
dir=`pwd`
echo $dir

rm final_results_RECALCULATED/BEDCOMPARE_SUMMARY_run_${run_ID}_N2.txt
rm final_results_RECALCULATED/FAMILY_TFPNs_F
for i in {1..30}
do
	i=$( echo "scale=4;$i / 10 " | bc -l )
	#FEED JUST SIMULATED to TEMP
	cd run_${run_ID}_filter_results_temp/
	cp run_${run_ID}_temp_nonredundant.bed SAVE.bed
	cat run_${run_ID}_temp_nonredundant.bed | awk -v i=$i '$7>=i {print $0}' > tmp && mv tmp run_${run_ID}_temp_nonredundant.bed
	python /lscr2/andersenlab/kml436/git_repos2/bamsurgeon/bed_compare9_redo_EDIT_RECALCULATE_absence.py 20000000 ${dir}/run_${run_ID}_filter_results_temp/ $TE_lengths adjusted_pos_SIM_tes_${run_ID}.bed run_${run_ID} $consensus_renamed >& run_${run_ID}_bedcompare.log
	mv SAVE.bed run_${run_ID}_temp_nonredundant.bed 
	cd ..


	####python /lscr2/andersenlab/kml436/git_repos2/bamsurgeon/bed_compare9_redo_EDIT_RECALCULATE_absence.py 20 /lscr2/andersenlab/kml436/RSV_TEST4/run_7/run_7_filter_results_temp/ /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt adjusted_pos_SIM_tes_7.bed run_7 /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/consensus.fasta >& run_7_bedcompare.log

	#FEED JUST REFERENCE to TELOCATE
	cp adjusted_pos_all_tes_${run_ID}.bed run_${run_ID}_filter_results_telocate/
	cd run_${run_ID}_filter_results_telocate/
	cp run_${run_ID}_telocate_nonredundant.bed SAVE.bed
	cat run_${run_ID}_telocate_nonredundant.bed | awk -v i=$i '$5>=i {print $0}' > tmp && mv tmp run_${run_ID}_telocate_nonredundant.bed
	python /lscr2/andersenlab/kml436/git_repos2/bamsurgeon/bed_compare9_redo_EDIT_RECALCULATE_absence.py 20000000 ${dir}/run_${run_ID}_filter_results_telocate/ $TE_lengths adjusted_pos_REF_tes_${run_ID}.bed run_${run_ID} $consensus_renamed >& run_${run_ID}_bedcompare.log
	mv SAVE.bed run_${run_ID}_telocate_nonredundant.bed
	cd ..


	#python /lscr2/andersenlab/kml436/git_repos2/bamsurgeon/bed_compare9_redo_EDIT_RECALCULATE_absence.py 20 /lscr2/andersenlab/kml436/RSV_TEST4/run_7/run_7_filter_results_telocate/ /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/fake_lengths.txt adjusted_pos_REF_tes_7.bed run_7 /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 


	#process TEMP results
	mkdir final_results_RECALCULATED
	cd run_${run_ID}_filter_results_temp/
	mkdir BEDCOMPARE_files_RECALCULATED
	cd intermediates_RECALCULATED/
	cat FAMILY_TFPN* > FAMILY_TFPNs_F_temp
	mv FAMILY_TFPNs_F_temp ../../final_results_RECALCULATED
	cd ..
	mv BEDCOMPARE_run_${run_ID}.txt BEDCOMPARE_FAMILY_run_${run_ID}.txt BEDCOMPARE_SUMMARY_run_${run_ID}.txt BEDCOMPARE_files_RECALCULATED
	cd BEDCOMPARE_files_RECALCULATED
	cat BEDCOMPARE_SUMMARY_run_${run_ID}.txt | sed 1d | sed 's/redundant\.bed/redundant\.bed_F/g' > tmp_filter_file.bed
	cd ../..

	#process TELOCATE results
	cd run_${run_ID}_filter_results_telocate/
	mkdir BEDCOMPARE_files_RECALCULATED
	cd intermediates_RECALCULATED/
	cat FAMILY_TFPN* > FAMILY_TFPNs_F_telocate
	mv FAMILY_TFPNs_F_telocate ../../final_results_RECALCULATED
	cd ..
	mv BEDCOMPARE_run_${run_ID}.txt BEDCOMPARE_FAMILY_run_${run_ID}.txt BEDCOMPARE_SUMMARY_run_${run_ID}.txt BEDCOMPARE_files_RECALCULATED
	cd BEDCOMPARE_files_RECALCULATED
	cat BEDCOMPARE_SUMMARY_run_${run_ID}.txt | sed 1d | sed 's/redundant\.bed/redundant\.bed_F/g' > tmp_filter_file.bed
	cd ../..


	#merge results
	max_dist=`cat run_${run_ID}_filter_results_temp/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed| cut -f3| sort -n| tail -n1`
	cat run_${run_ID}_filter_results_temp/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed | awk  -v i=$i -v max_dist=$max_dist '$3==max_dist {print $0"\t"i}' > tmp && mv tmp run_${run_ID}_filter_results_temp/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed
	max_dist=`cat run_${run_ID}_filter_results_telocate/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed| cut -f3| sort -n| tail -n1`
	cat run_${run_ID}_filter_results_telocate/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed | awk  -v i=$i -v max_dist=$max_dist '$3==max_dist {print $0"\t"i}' > tmp && mv tmp run_${run_ID}_filter_results_telocate/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed
	cat run_${run_ID}_filter_results_telocate/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed | awk  -v i=$i '{print $0"\t"i}' > tmp && mv tmp run_${run_ID}_filter_results_telocate/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed 
	cat run_${run_ID}_filter_results_temp/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed run_${run_ID}_filter_results_telocate/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed > BEDCOMPARE_SUMMARY_run_${run_ID}_N2.txt
	cat BEDCOMPARE_SUMMARY_run_${run_ID}_N2.txt >> final_results_RECALCULATED/BEDCOMPARE_SUMMARY_run_${run_ID}_N2.txt
	cd final_results_RECALCULATED
	cat FAMILY_TFPNs_F_temp FAMILY_TFPNs_F_telocate | awk  -v i=$i '{print $0"\t"i}' >> FAMILY_TFPNs_F

	#flip distnace and read support
	#remove summary files before hand
	cd ..
done
cd final_results_RECALCULATED
#flip distance and read support columns to put in compatible format for mean-calculating script
cat BEDCOMPARE_SUMMARY_run_${run_ID}_N2.txt | awk '{print $1"\t"$2"\t"$7"\t"$4"\t"$5"\t"$6"\t"$3}' > tmp && mv tmp BEDCOMPARE_SUMMARY_run_${run_ID}_N2.txt
echo "Done"