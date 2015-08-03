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
bam=/lscr2/andersenlab/dec211/RUN/v2_snpset/bam/JU1395.bam
TE_consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
processors=18
TTR=/lscr2/andersenlab/kml436/git_repos2/mcclintock/
#HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_tes.gff ###
HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names_alias.bed
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
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/AB-PR/consensus_wTC8.fasta 
bam_name=$(basename "$bam" .bam)
echo "bam name is ${bam_name}"
mkdir ${run_ID}_${bam_name}
dir1=`pwd`
cd ${run_ID}_${bam_name}
dir=`pwd`
##############################################################
##create error log file for everything?
##randomsites--add options here?


echo "Generating random intervals to attempt spike ins......"
$python_version ${bam_surgeon}/etc/randomsites.py -g ${reference} -n 10 --minvaf 1 --maxvaf 1 sv > ${run_ID}_BS_random_intervals.txt
bedtools intersect -v -a ${run_ID}_BS_random_intervals.txt -b $high_coverage_regions > interval_tmp && mv interval_tmp ${run_ID}_BS_random_intervals.txt
# avoid spiking in to already existing transposons
bedtools intersect -v -a ${run_ID}_BS_random_intervals.txt -b $existing_transposons > interval_tmp && mv interval_tmp ${run_ID}_BS_random_intervals.txt

##spike ins (using altered csripts)y
echo "Spiking in mutations......"
$python_version ${bam_surgeon}/addsv.py -v ${run_ID}_BS_random_intervals.txt -f $bam -r $reference -o ${run_ID}_${bam_name}.bam --inslib $all_WB_seqs -p $processors -l 600 >& ${run_ID}_${bam_name}_bamsurgeon.log
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
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' "${run_ID}_${bam_name}_temp_sorted_precut_raw.txt" >> "${run_ID}_${bam_name}_temp_raw.bed"
mv ${run_ID}_${bam_name}_temp_raw.bed ${run_ID}_${bam_name}_non_filter_results/
echo "Filtering TEMP output......"
# Filter for insertions with support for both ends
awk '{if ($8 == "1p1") print $0}' "${run_ID}_${bam_name}_temp_sorted_precut_raw.txt" > "${run_ID}_${bam_name}_temp_sorted_redundant.txt"
cut -f1-3,5-7 "${run_ID}_${bam_name}_temp_sorted_redundant.txt" >> "${run_ID}_${bam_name}_temp_redundant.bed"
# Filter out redundant insertions
sort -k1,3 -k4rn "${run_ID}_${bam_name}_temp_sorted_redundant.txt" | sort -u -k1,3 | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' > "${run_ID}_${bam_name}_tmp"
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
echo "running RETROSEQ-discovery..."
##Discovery
perl ${TTR}/RetroSeq/bin/retroseq.pl -discover -id 70 -bam ${run_ID}_${bam_name}.sorted.bam -eref $element_list -align -output ${run_ID}_${bam_name}_align.discovery
echo "running RETROSEQ-calling..."
##Calling
perl ${TTR}/RetroSeq/bin/retroseq.pl -call -bam ${run_ID}_${bam_name}.sorted.bam -input ${run_ID}_${bam_name}_align.discovery -filter $location_list -ref $reference -output ${run_ID}_${bam_name}_align.calling -orientate yes
echo "RETROSEQ finished...altering output..."
awk '$1!~/#/{print $0}' "${run_ID}_${bam_name}_align.calling.PE.vcf" > "Retro/tmp"
awk -F'[=,\t:]' '{print $1"\t"$11"\t"$12"\t"$22"\t"$10"_non-reference_retroseq_rp_\t0\t.\t"$21}' "Retro/tmp" > "${run_ID}_${bam_name}_align_retroseq_presort.txt"
bedtools sort -i "${run_ID}_${bam_name}_align_retroseq_presort.txt" > "${run_ID}_${bam_name}_align_retroseq_sorted_raw.txt"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7"\t"$8}' "${run_ID}_${bam_name}_align_retroseq_sorted_raw.txt" > "${run_ID}_${bam_name}_align_retroseq_sorted_counted_raw.txt"
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' "${run_ID}_${bam_name}_align_retroseq_sorted_counted_raw.txt" >> "${run_ID}_${bam_name}_align_retroseq_raw.bed"
mv ${run_ID}_${bam_name}_align_retroseq_raw.bed ${run_ID}_${bam_name}_non_filter_results/
## Filter results
echo "Filtering RETROSEQ output......"
awk '{if ($8 >= 6) print $0}' "${run_ID}_${bam_name}_align_retroseq_sorted_counted_raw.txt" > "${run_ID}_${bam_name}_align_retroseq_redundant.txt"
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' "${run_ID}_${bam_name}_align_retroseq_redundant.txt" >> "${run_ID}_${bam_name}_align_retroseq_redundant.bed"
## Filter for redundant predictions
sort -k1,3 -k4rn "${run_ID}_${bam_name}_align_retroseq_redundant.txt" | sort -u -k1,3 | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' > Retro/temp
bedtools sort -i Retro/temp >> "${run_ID}_${bam_name}_align_retroseq_nonredundant.bed"
mv ${run_ID}_${bam_name}_align_retroseq_nonredundant.bed ${run_ID}_${bam_name}_filter_results/

cd ${run_ID}_${bam_name}_non_filter_results
mv ${run_ID}_${bam_name}_align_retroseq_raw.bed ${run_ID}_${bam_name}_retroseq_nonredundant.bed
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_retroseq_nonredundant.bed TMP
mv ${run_ID}_${bam_name}_retroseq_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy.bed
mv TMP ${run_ID}_${bam_name}_retroseq_nonredundant.bed
cat ${run_ID}_${bam_name}_retroseq_nonredundant.bed| sed s'/chrX/X/g'| sed '/^track/d' | awk '$4 ~/non-reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${bam_name}_temp && mv ${run_ID}_${bam_name}_retroseq_nonredundant.bed 

cd ..
cd ${run_ID}_${bam_name}_filter_results
mv ${run_ID}_${bam_name}_align_retroseq_nonredundant.bed ${run_ID}_${bam_name}_retroseq_nonredundant.bed
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_retroseq_nonredundant.bed TMP
mv ${run_ID}_${bam_name}_retroseq_nonredundant.bed ${run_ID}_${bam_name}_nonredundant_org_copy.bed
mv TMP ${run_ID}_${bam_name}_retroseq_nonredundant.bed
cat ${run_ID}_${bam_name}_retroseq_nonredundant.bed| sed s'/chrX/X/g'| sed '/^track/d' | awk '$4 ~/non-reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${bam_name}_temp && mv ${run_ID}_${bam_name}_retroseq_nonredundant.bed 
cd ..

#########################
echo "Creating sam file..."
mkdir ${run_ID}_${bam_name}_sam_file
samtools view ${run_ID}_${bam_name}.sorted.bam |sort --temporary-directory=${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_sam_file > ${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_sam_file/${run_ID}_${bam_name}.sorted.sam
##in  /lscr2/andersenlab/kml436/git_repos2/mcclintock/TE-locate/
echo "running TELOCATE..."

mkdir TELOCATE
cd ${TTR}/TE-locate/
perl TE_locate.pl 2 ${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_sam_file/ $HL_gff $reference ${dir1}/${run_ID}_${bam_name}/TELOCATE/TEL $minimal_Distance_to_count $minimal_supporting_reads $minimal_supporting_individuals
cd ${dir1}/${run_ID}_${bam_name}
echo "TELOCATE finished...altering output..."
##LEFT OFF HERE
sed -e '1,2d' "${dir1}/${run_ID}_${bam_name}/TELOCATE/TEL_${minimal_Distance_to_count}_reads${minimal_supporting_reads}_acc${minimal_supporting_individuals}.info" > "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile"
awk -F'[\t/]' '{printf $1"\t"; if($17=="old") printf $2-1"\t"$2"\t"$8"\t"$5"_reference_telocate_rp_\t0"; else if($17=="new") printf $2-1"\t"$2"\t"$8"\t"$5"_non-reference_telocate_rp_\t0"; if($14=="parallel") print "\t+"; else if($14 =="uncertain") print "\t."; else print "\t-";}' "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile" > "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_presort.txt"
bedtools sort -i "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_presort.txt" > "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_sorted_redundant.txt"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7}' "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_sorted_redundant.txt" > "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_counted_redundant.txt"
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_counted_redundant.txt" >> "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed"
cat ${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed | awk '$4~/_non-reference/ {print $0}' > tmp && mv tmp ${run_ID}_tmpfile_telocate_redundant.bed 
mv ${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed ${run_ID}_${bam_name}_non_filter_results/
echo "STEP1"
cd ${run_ID}_${bam_name}_non_filter_results/
#rename TELOCATE NON FILTER results (from WB element names to consensus names)....won't actualy end up using the "unfiltered" results...same as "filtered"
mv ${run_ID}_${bam_name}_tmpfile_telocate_redundant.bed ${run_ID}_${bam_name}_telocate_nonredundant.bed 
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_telocate_nonredundant.bed  TMP
mv ${run_ID}_${bam_name}_telocate_nonredundant.bed  ${run_ID}_nonredundant_org_copy.bed
mv TMP ${run_ID}_${bam_name}_telocate_nonredundant.bed 
!!!
closestBed -a !$ref_tes! -b ${run_ID}_${bam_name}_telocate_nonredundant.bed   -d -t all | awk '$13<=1000 {print $0}' > ${bam_name}_closest.txt
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/process_telocate_output.py ${bam_name}_closest.txt
cat processed_calls.txt| sort -k1,1 -k2,2n > ${run_ID}_${bam_name}_telocate_nonredundant.bed 
!!!
cat ${run_ID}_${bam_name}_telocate_nonredundant.bed| awk '$4 ~/_non-reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${run_ID}_temp && mv ${run_ID}_temp ${run_ID}_${bam_name}_telocate_nonredundant.bed 
cd ..

echo "Filtering TELOCATE output......"
sort -k1,3 -k4rn "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfile_telocate_counted_redundant.txt"| sort -u -k1,3 | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' > "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfileNR"
bedtools sort -i "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_tmpfileNR" >> "${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed"
cat ${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed| awk '$4~/_non-reference/ {print $0}' > tmp && mv tmp ${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed
mv ${dir1}/${run_ID}_${bam_name}/${run_ID}_${bam_name}_telocate_nonredundant.bed ${run_ID}_${bam_name}_filter_results/
######

#rename TELOCATE results (from WB element names to consensus names)
cd ${run_ID}_${bam_name}_filter_results/
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${run_ID}_${bam_name}_telocate_nonredundant.bed  TMP
mv ${run_ID}_${bam_name}_telocate_nonredundant.bed  ${run_ID}_nonredundant_org_copy.bed
mv TMP ${run_ID}_${bam_name}_telocate_nonredundant.bed
!!!
!!!
cat ${run_ID}_${bam_name}_telocate_nonredundant.bed| awk '$4 ~/_non-reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${run_ID}_temp && mv ${run_ID}_temp ${run_ID}_${bam_name}_telocate_nonredundant.bed 

cd ..
###################
##
#
#
#
#


echo "Comparing outputs to simulated structural variants......"
python /exports/people/andersenlab/kml436/scripts/create_bam_surgeon_log.py ${run_ID}_${bam_name}_bamsurgeon.log
bedtools sort -i ${run_ID}_${bam_name}_final_positions.bed> tmpSUM && mv tmpSUM ${run_ID}_${bam_name}_final_positions.bed
cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_final_positions_WB_copy.bed 
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2.py $family_renames ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_final_positions_TMP.bed
mv ${run_ID}_${bam_name}_final_positions_TMP.bed ${run_ID}_${bam_name}_final_positions.bed
cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_non_filter_results/
cp ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name}_filter_results/

mkdir final_results

cd ${run_ID}_${bam_name}_non_filter_results/
mkdir BEDCOMPARE_files
${python_version} ${bam_surgeon}/bed_compare9_redo_EDIT.py 20 ${dir}/${run_ID}_${bam_name}_non_filter_results/ $fake_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} $consensus_renamed >& ${run_ID}_${bam_name}_bedcompare.log
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
${python_version} ${bam_surgeon}/bed_compare9_redo_EDIT.py 20 ${dir}/${run_ID}_${bam_name}_filter_results/ $fake_lengths ${run_ID}_${bam_name}_final_positions.bed ${run_ID}_${bam_name} $consensus_renamed >& ${run_ID}_${bam_name}_bedcompare.log
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
