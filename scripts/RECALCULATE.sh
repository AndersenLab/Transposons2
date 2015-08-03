#!/bin/bash
#this script recalulcates the distances based on the starting position of the found transposons rather then the whole transposons (unless the TE
# was found as the exact position of the simualted TE , the distance will not be 0)
# (runs bed_compare9_redo_RECALCULATE)
# USE: RECALCUALTE.sh
# or USE: RECALCUALTE.sh <length_file>

bam=/lscr2/andersenlab/dec211/RUN/v2_snpset/bam/N2.bam
bam_name=$(basename "$bam" .bam)
echo "bam name is ${bam_name}"
python_version=/exports/people/andersenlab/kml436/python2/bin/python
bam_surgeon=/lscr2/andersenlab/kml436/git_repos2/bamsurgeon
TE_lengths=${1}
#TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/round3_lengths.txt
#TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/lengths.txt
#TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/WB_all_lengths/lengths.txt


for i in {1..8}
do
	run_ID=run_${i}
	cd ${run_ID}_${bam_name}
	dir=`pwd`
	cd ${run_ID}_${bam_name}_non_filter_results/
	mkdir BEDCOMPARE_files_RECALCULATED
	${python_version} ${bam_surgeon}/bed_compare9_redo_RECALCULATE.py 20 ${dir}/${run_ID}_${bam_name}_non_filter_results/ $TE_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} >& ${run_ID}_${bam_name}_bedcompare.log
	mv BEDCOMPARE_${run_ID}_${bam_name}.txt BEDCOMPARE_FAMILY_${run_ID}_${bam_name}.txt BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt BEDCOMPARE_files_RECALCULATED
	cd BEDCOMPARE_files_RECALCULATED
	cat BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt | sed 's/redundant\.bed/redundant\.bed_NF/g' > tmp_non_filter_file.bed
	cd ../..

	cd ${run_ID}_${bam_name}_filter_results/
	mkdir BEDCOMPARE_files_RECALCULATED
	${python_version} ${bam_surgeon}/bed_compare9_redo_RECALCULATE.py 20 ${dir}/${run_ID}_${bam_name}_non_filter_results/ $TE_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} >& ${run_ID}_${bam_name}_bedcompare.log
	mv BEDCOMPARE_${run_ID}_${bam_name}.txt BEDCOMPARE_FAMILY_${run_ID}_${bam_name}.txt BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt BEDCOMPARE_files_RECALCULATED
	cd BEDCOMPARE_files_RECALCULATED
	cat BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt | sed 's/redundant\.bed/redundant\.bed_F/g' > tmp_filter_file.bed
	cd ../..

	cat ${run_ID}_${bam_name}_non_filter_results/BEDCOMPARE_files_RECALCULATED/tmp_non_filter_file.bed ${run_ID}_${bam_name}_filter_results/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed > BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt
	mkdir final_results_RECALCULATED
	mv BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt final_results_RECALCULATED
	cd ..

done



