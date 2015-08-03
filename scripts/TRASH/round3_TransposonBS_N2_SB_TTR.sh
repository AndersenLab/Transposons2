#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=18000

python_version=/exports/people/andersenlab/kml436/python2/bin/python
bam_surgeon=/lscr2/andersenlab/kml436/git_repos2/bamsurgeon
reference=/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa
run_ID=run_${1}
bam=/lscr2/andersenlab/dec211/RUN/v2_snpset/bam/N2.bam
TE_consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round3_consensus_fasta.fasta
##ALTERNATE: TE_consensus=/lscr2/andersenlab/kml436/common.fasta
processors=18
TTR=/lscr2/andersenlab/kml436/git_repos2/mcclintock/
element_list=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round3_retroseq_files/fasta_list.txt
##ALTERNATE: element_list=/lscr2/andersenlab/kml436/sv_files/WB_REPB_elementlist.gff
location_list=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round3_retroseq_files/bed_list.txt
##ALTERNATE: location_list=/lscr2/andersenlab/kml436/sv_files/WB_REPB_locationlist.gff
HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round2_WB_tes.gff #same for round3
##ALTERNATE: HL_gff=/lscr2/andersenlab/kml436/sv_files/WB_REPB_ID_transposable_element_HL.gff
minimal_Distance_to_count=1000
minimal_supporting_reads=3
minimal_supporting_individuals=1
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round3_lengths.txt
high_coverage_regions=/lscr2/andersenlab/kml436/sv_files/N2_high_coverage_regions
existing_transposons=/lscr2/andersenlab/kml436/ID_transposable_element.bed
bam_name=$(basename "$bam" .bam)
echo "bam name is ${bam_name}"
mkdir ${run_ID}_${bam_name}
cd ${run_ID}_${bam_name}
dir=`pwd`
##############################################################
##create error log file for everything?
##randomsites--add options here?


echo "Generating random intervals to attempt spike ins......"
$python_version ${bam_surgeon}/etc/randomsites.py -g ${reference} -n 1000 --minvaf 1 --maxvaf 1 sv > ${run_ID}_BS_random_intervals.txt
bedtools intersect -v -a ${run_ID}_BS_random_intervals.txt -b $high_coverage_regions > interval_tmp && mv interval_tmp ${run_ID}_BS_random_intervals.txt
# avoid spiking in to already existing transposons
bedtools intersect -v -a ${run_ID}_BS_random_intervals.txt -b $existing_transposons > interval_tmp && mv interval_tmp ${run_ID}_BS_random_intervals.txt

##spike ins (using altered csripts)y
echo "Spiking in mutations......"
$python_version ${bam_surgeon}/addsv.py -v ${run_ID}_BS_random_intervals.txt -f $bam -r $reference -o ${run_ID}_${bam_name}.bam --inslib $TE_consensus -p $processors -l 600 >& ${run_ID}_${bam_name}_bamsurgeon.log
##sort bam files
echo "Sorting simulated bam file..."
/opt/samtools/x64/samtools/samtools sort -o -@ 8 ${run_ID}_${bam_name}.bam out > ${run_ID}_${bam_name}.sorted.bam
/opt/samtools/x64/samtools/samtools index ${run_ID}_${bam_name}.sorted.bam
##create sam files

mkdir ${run_ID}_${bam_name}_non_filter_results/
mkdir ${run_ID}_${bam_name}_filter_results/
mkdir Retro
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
echo "running RETROSEQ-discovery..."
##Discovery
perl ${TTR}/RetroSeq/bin/retroseq.pl -discover -bam ${run_ID}_${bam_name}.sorted.bam -eref $element_list -align -output ${run_ID}_${bam_name}_align.discovery
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
##############################################################
echo "Creating sam file..."
mkdir ${run_ID}_${bam_name}_sam_file
samtools view ${run_ID}_${bam_name}.sorted.bam |sort --temporary-directory=${dir}/${run_ID}_${bam_name}_sam_file > ${dir}/${run_ID}_${bam_name}_sam_file/${run_ID}_${bam_name}.sorted.sam
#in  /lscr2/andersenlab/kml436/git_repos2/mcclintock/TE-locate/
echo "running TELOCATE..."
mkdir ${dir}/${run_ID}_${bam_name}
mkdir ${dir}/${run_ID}_${bam_name}/TELOCATE
cd ${TTR}/TE-locate/
perl TE_locate.pl 2 ${dir}/${run_ID}_${bam_name}_sam_file/ $HL_gff $reference ${dir}/${run_ID}_${bam_name}/TELOCATE/TEL $minimal_Distance_to_count $minimal_supporting_reads $minimal_supporting_individuals
cd ${dir}
echo "TELOCATE finished...altering output..."
##LEFT OFF HERE
sed -e '1,2d' "${dir}/${run_ID}_${bam_name}/TELOCATE/TEL_${minimal_Distance_to_count}_reads${minimal_supporting_reads}_acc${minimal_supporting_individuals}.info" > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile"
awk -F'[\t/]' '{printf $1"\t"; if($17=="old") printf $2-1"\t"$2"\t"$8"\t"$5"_reference_telocate_rp_\t0"; else if($17=="new") printf $2-1"\t"$2"\t"$8"\t"$5"_non-reference_telocate_rp_\t0"; if($14=="parallel") print "\t+"; else if($14 =="uncertain") print "\t."; else print "\t-";}' "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile" > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_presort.txt"
bedtools sort -i "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_presort.txt" > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_sorted_redundant.txt"

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7}' "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_sorted_redundant.txt" > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_counted_redundant.txt"
cut -f1-3,5-7 "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_counted_redundant.txt" >> "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed"
mv ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed ${run_ID}_${bam_name}_non_filter_results/
echo "Filtering TELOCATE output......"
sort -k1,3 -k4rn "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_counted_redundant.txt"| sort -u -k1,3 | cut -f1-3,5-7 > "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfileNR"
bedtools sort -i "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfileNR" >> "${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed"
mv ${dir}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed ${run_ID}_${bam_name}_filter_results/
######
cd ${run_ID}_${bam_name}_non_filter_results
mv ${run_ID}_${bam_name}_temp_raw.bed  ${run_ID}_${bam_name}_temp_nonredundant.bed
mv ${run_ID}_${bam_name}_align_retroseq_raw.bed ${run_ID}_${bam_name}_retroseq_nonredundant.bed
mv ${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed ${run_ID}_${bam_name}_telocate_nonredundant.bed
cd ..
cd ${run_ID}_${bam_name}_filter_results
mv ${run_ID}_${bam_name}_align_retroseq_nonredundant.bed ${run_ID}_${bam_name}_retroseq_nonredundant.bed
cd ..
#########################
echo "Comparing outputs to simulated structural variants......"
python /exports/people/andersenlab/kml436/scripts/create_bam_surgeon_log.py ${run_ID}_${bam_name}_bamsurgeon.log
bedtools sort -i ${run_ID}_${bam_name}_final_positions.bed> tmpSUM && mv tmpSUM ${run_ID}_${bam_name}_final_positions.bed
cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_non_filter_results/
cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_filter_results/

cd ${run_ID}_${bam_name}_non_filter_results/
mkdir BEDCOMPARE_files
${python_version} ${bam_surgeon}/bed_compare9_redo.py 100 ${dir}/${run_ID}_${bam_name}_non_filter_results/ $TE_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} >& ${run_ID}_${bam_name}_bedcompare.log
mv BEDCOMPARE_${run_ID}_${bam_name}.txt BEDCOMPARE_FAMILY_${run_ID}_${bam_name}.txt BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt BEDCOMPARE_files
cd BEDCOMPARE_files
cat BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt | sed 's/redundant\.bed/redundant\.bed_NF/g' > tmp_non_filter_file.bed
cd ../..

cd ${run_ID}_${bam_name}_filter_results/
mkdir BEDCOMPARE_files
${python_version} ${bam_surgeon}/bed_compare9_redo.py 100 ${dir}/${run_ID}_${bam_name}_filter_results/ $TE_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} >& ${run_ID}_${bam_name}_bedcompare.log
mv BEDCOMPARE_${run_ID}_${bam_name}.txt BEDCOMPARE_FAMILY_${run_ID}_${bam_name}.txt BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt BEDCOMPARE_files
cd BEDCOMPARE_files
cat BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt | sed 1d | sed 's/redundant\.bed/redundant\.bed_F/g' > tmp_filter_file.bed
cd ../..

cat ${run_ID}_${bam_name}_non_filter_results/BEDCOMPARE_files/tmp_non_filter_file.bed ${run_ID}_${bam_name}_filter_results/BEDCOMPARE_files/tmp_filter_file.bed > BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt
mkdir final_results
mv BEDCOMPARE_SUMMARY_${run_ID}_${bam_name}.txt final_results

#TEMP(){
#	echo "temp function is $1"
#}
#bash ${TTR}/TEMP/scripts/TEMP_Insertion.sh -i ${run_ID}_${bam_name}.sorted.bam -s ${TTR}/TEMP/scripts/ -x 30 -r $TE_consensus -m 1 -c 20 &> ${run_ID}_TEMP_log.txt
#samtools version
#bedcompare
##########
