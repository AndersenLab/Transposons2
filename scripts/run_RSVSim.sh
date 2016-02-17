#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --exclude=node[6,8]
# this script runs the RSVSimtransposon simulation and detection scripts 
# USE: run_RSVSim.sh <run_number>
consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta
#consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/consensus.fasta
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths_plusTC8.txt
#TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt
WB_elements=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/Wb-TC8.fasta #these are the plus strands without TC8
family_renames=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt

##Rockman was not included in the N2 bam, so do not include it here
##A_N2_1=/lscr2/andersenlab/kml436/fasta_subset/Rockman096-R096.2-N2-paired-1.fq
##A_N2_2=/lscr2/andersenlab/kml436/fasta_subset/Rockman096-R096.2-N2-paired-2.fq
#B_N2_1=/lscr2/andersenlab/kml436/fasta_subset/Princeton-P-N2_CGC-paired-1.fq
#B_N2_2=/lscr2/andersenlab/kml436/fasta_subset/Princeton-P-N2_CGC-paired-2.fq
#C_N2_1=/lscr2/andersenlab/kml436/fasta_subset/Uchicago-L001-N2Baer-paired-1.fq
#C_N2_2=/lscr2/andersenlab/kml436/fasta_subset/Uchicago-L001-N2Baer-paired-2.fq

###
#D_N2_1=/lscr2/andersenlab/kml436/fasta_subset/BGI1-RET1-N2_CGC-paired-1.fq
#D_N2_2=/lscr2/andersenlab/kml436/fasta_subset/BGI1-RET1-N2_CGC-paired-2.fq
#E_N2_1=/lscr2/andersenlab/kml436/fasta_subset/BGI1-RET6-N2_HRH-paired-1.fq
#E_N2_2=/lscr2/andersenlab/kml436/fasta_subset/BGI1-RET6-N2_HRH-paired-2.fq
#F_N2_1=/lscr2/andersenlab/kml436/fasta_subset/BGI2-RET1-N2_CGC-paired-1.fq
#F_N2_2=/lscr2/andersenlab/kml436/fasta_subset/BGI2-RET1-N2_CGC-paired-2.fq
#G_N2_1=/lscr2/andersenlab/kml436/fasta_subset/BGI2-RET6-N2_HRH-paired-1.fq
#G_N2_2=/lscr2/andersenlab/kml436/fasta_subset/BGI2-RET6-N2_HRH-paired-2.fq
#H_N2_1=/lscr2/andersenlab/kml436/fasta_subset/BGI3-RET1a-N2_CGC-paired-1.fq
#H_N2_2=/lscr2/andersenlab/kml436/fasta_subset/BGI3-RET1a-N2_CGC-paired-2.fq
#I_N2_1=/lscr2/andersenlab/kml436/fasta_subset/BGI3-RET6b-N2_HRH-paired-1.fq
#I_N2_2=/lscr2/andersenlab/kml436/fasta_subset/BGI3-RET6b-N2_HRH-paired-2.fq


TEMP_scripts=/lscr2/andersenlab/kml436/git_repos2/mcclintock/TEMP/scripts/
absence=/lscr2/andersenlab/kml436/git_repos2/mcclintock/TEMP/scripts/TEMP_Absence.sh
faToTwoBit_dir=/exports/people/andersenlab/kml436
original_ref_pos=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names.gff
#original_ref_pos=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_consensus_RSV_sim.bed
TTR=/lscr2/andersenlab/kml436/git_repos2/mcclintock/
minimal_Distance_to_count=1000
minimal_supporting_reads=1
minimal_supporting_individuals=1
run_ID=${1}
mkdir run_${run_ID}
cd run_${run_ID}
dir=`pwd`
echo $dir

find . -name "*unpair*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*insertion*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*align*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*uniq*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*clipped*" -type f -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*filter_results*" -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;
find . -name "*final_results*" -exec rename 's/(\w+)$/$1\_IGNORE/' {} \;

mkdir run_${run_ID}_filter_results_temp/
mkdir run_${run_ID}_filter_results_telocate/




##choose 100 TEs from the consensus TE seqs to randomly simulate in the reference
#python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/create_RSVSim_genome.py  $WB_elements $run_ID

##run RSV SIM
#Rscript /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/RSV_SIM_OLDER6.R RSVSIM_genome_${run_ID}.fasta RSVSIM_te_pos_${run_ID}.txt
#mv insertions.csv insertions_${run_ID}.csv
#cp insertions_${run_ID}.csv original_insertion_file_${run_ID}.txt
#cat insertions_${run_ID}.csv | sed 's/_copy[1-9]//g' | awk '{print $5"\t"$6"\t"$7"\t"$2"\t0\t+"}'|sort -k1,1 -k2,2n > tmp && mv tmp insertions_${run_ID}.csv

##from the RSV genome, pull out only the original chromosomes containing the simualted TEs (removes the "copy" TEs taken in addition to the original genome)
#python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/edit_rearranged_genome.py genome_rearranged.fasta ${run_ID}
#sed '1d' insertions_${run_ID}.csv > insertions_${run_ID}.bed
#rm insertions_${run_ID}.csv
#cat $original_ref_pos insertions_${run_ID}.bed > all_tes_${run_ID}.bed
#cat all_tes_${run_ID}.bed|sort -k1,1 -k2,2n > tmp && mv tmp all_tes_${run_ID}.bed

##adjust the postions of the TEs to account for the simulated TEs
#python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/adjust_RSVSIM_positions.py all_tes_${run_ID}.bed

##rename the position file (from WB element names to consensus names)
#cp adjusted_pos_all_tes_${run_ID}.bed  adjusted_pos_all_tes_${run_ID}_WB_copy.bed
#python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2.py $family_renames adjusted_pos_all_tes_${run_ID}.bed adjusted_pos_all_tes_${run_ID}_TMP.bed
#mv adjusted_pos_all_tes_${run_ID}_TMP.bed adjusted_pos_all_tes_${run_ID}.bed

##created position file usable as input to telocate
cat adjusted_pos_all_tes_${run_ID}.bed | awk '{print $1"\tTransposon\t"$4"\t"$2+1"\t"$3+1"\t"$5"\t"$6"\t.\tID="$4";Name="$4";Alias="$4}' > telocate_adjusted_pos_all_tes_${run_ID}.bed

#cat adjusted_pos_all_tes_${run_ID}.bed | awk '{print $1"\tTransposon\t"$4"\t"$2"\t"$3"\t"$5"\t"$6"\t.\tID="$4";Name="$4";Alias="$4}' > telocate_adjusted_pos_all_tes_${run_ID}.bed
#/opt/samtools/x64/samtools/samtools faidx new_genome_${run_ID}.fasta
#bwa index new_genome_${run_ID}.fasta

##BWA MEM
##bwa mem -t 24 -R '@RG\tID:A_N2\tSM:N2' new_genome_${run_ID}.fasta $A_N2_1 $A_N2_2 > A_N2.sam
#bwa mem -t 24 -R '@RG\tID:B_N2\tSM:N2' new_genome_${run_ID}.fasta $B_N2_1 $B_N2_2 -M > B_N2.sam
#bwa mem -t 24 -R '@RG\tID:C_N2\tSM:N2' new_genome_${run_ID}.fasta $C_N2_1 $C_N2_2 -M > C_N2.sam
#bwa mem -t 24 -R '@RG\tID:D_N2\tSM:N2' new_genome_${run_ID}.fasta $D_N2_1 $D_N2_2 -M > D_N2.sam
#bwa mem -t 24 -R '@RG\tID:E_N2\tSM:N2' new_genome_${run_ID}.fasta $E_N2_1 $E_N2_2 -M > E_N2.sam
#bwa mem -t 24 -R '@RG\tID:F_N2\tSM:N2' new_genome_${run_ID}.fasta $F_N2_1 $F_N2_2 -M > F_N2.sam
#bwa mem -t 24 -R '@RG\tID:G_N2\tSM:N2' new_genome_${run_ID}.fasta $G_N2_1 $G_N2_2 -M > G_N2.sam
#bwa mem -t 24 -R '@RG\tID:H_N2\tSM:N2' new_genome_${run_ID}.fasta $H_N2_1 $H_N2_2 -M > H_N2.sam
#bwa mem -t 24 -R '@RG\tID:I_N2\tSM:N2' new_genome_${run_ID}.fasta $I_N2_1 $I_N2_2 -M > I_N2.sam



#SAMTOOLS VIEW
#/opt/samtools/x64/samtools/samtools view A_N2.sam -@ 24 -bhu -o A_N2.unsorted.bam
#/opt/samtools/x64/samtools/samtools view B_N2.sam -@ 24 -bhu -o B_N2.unsorted.bam
#/opt/samtools/x64/samtools/samtools view C_N2.sam -@ 24 -bhu -o C_N2.unsorted.bam
#/opt/samtools/x64/samtools/samtools view D_N2.sam -@ 24 -bhu -o D_N2.unsorted.bam
#/opt/samtools/x64/samtools/samtools view E_N2.sam -@ 24 -bhu -o E_N2.unsorted.bam
#/opt/samtools/x64/samtools/samtools view F_N2.sam -@ 24 -bhu -o F_N2.unsorted.bam
#/opt/samtools/x64/samtools/samtools view G_N2.sam -@ 24 -bhu -o G_N2.unsorted.bam
#/opt/samtools/x64/samtools/samtools view H_N2.sam -@ 24 -bhu -o H_N2.unsorted.bam
#/opt/samtools/x64/samtools/samtools view I_N2.sam -@ 24 -bhu -o I_N2.unsorted.bam

#SAMTOOLS SORT
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp1 A_N2.unsorted.bam > A_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp2 B_N2.unsorted.bam > B_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp3 C_N2.unsorted.bam > C_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp4 D_N2.unsorted.bam > D_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp5 E_N2.unsorted.bam > E_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp6 F_N2.unsorted.bam > F_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp7 G_N2.unsorted.bam > G_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp8 H_N2.unsorted.bam > H_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp9 I_N2.unsorted.bam > I_N2.sorted.bam


#SAMTOOLS INDEX
#/opt/samtools/x64/samtools/samtools index A_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index B_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index C_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index D_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index E_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index F_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index G_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index H_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index I_N2.sorted.bam

##PICARD MARK DUPLICATES
#java -jar /opt/picard-tools/MarkDuplicates.jar I=B_N2.sorted.bam O=MD_B_N2.sorted.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=B_N2_metrics.txt ASSUME_SORTED=true
#mv MD_B_N2.sorted.bam B_N2.sorted.bam
#java -jar /opt/picard-tools/MarkDuplicates.jar I=C_N2.sorted.bam O=MD_C_N2.sorted.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=C_N2_metrics.txt ASSUME_SORTED=true
#mv MD_C_N2.sorted.bam C_N2.sorted.bam
#java -jar /opt/picard-tools/MarkDuplicates.jar I=D_N2.sorted.bam O=MD_D_N2.sorted.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=D_N2_metrics.txt ASSUME_SORTED=true
#mv MD_D_N2.sorted.bam D_N2.sorted.bam
#java -jar /opt/picard-tools/MarkDuplicates.jar I=E_N2.sorted.bam O=MD_E_N2.sorted.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=E_N2_metrics.txt ASSUME_SORTED=true
#mv MD_E_N2.sorted.bam E_N2.sorted.bam
#java -jar /opt/picard-tools/MarkDuplicates.jar I=F_N2.sorted.bam O=MD_F_N2.sorted.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=F_N2_metrics.txt ASSUME_SORTED=true
#mv MD_F_N2.sorted.bam F_N2.sorted.bam
#java -jar /opt/picard-tools/MarkDuplicates.jar I=G_N2.sorted.bam O=MD_G_N2.sorted.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=G_N2_metrics.txt ASSUME_SORTED=true
#mv MD_G_N2.sorted.bam G_N2.sorted.bam
#java -jar /opt/picard-tools/MarkDuplicates.jar I=H_N2.sorted.bam O=MD_H_N2.sorted.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=H_N2_metrics.txt ASSUME_SORTED=true
#mv MD_H_N2.sorted.bam H_N2.sorted.bam
#java -jar /opt/picard-tools/MarkDuplicates.jar I=I_N2.sorted.bam O=MD_I_N2.sorted.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=I_N2_metrics.txt ASSUME_SORTED=true
#mv MD_I_N2.sorted.bam I_N2.sorted.bam


##SAMTOOLS RE_INDEX (not necessary?)
##/opt/samtools/x64/samtools/samtools index A_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index B_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index C_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index D_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index E_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index F_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index G_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index H_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index I_N2.sorted.bam


##SAMTOOLS MERGE AND INDEX
#/opt/samtools/x64/samtools/samtools merge -f merged_bam_${run_ID}.sorted.bam B_N2.sorted.bam C_N2.sorted.bam D_N2.sorted.bam E_N2.sorted.bam F_N2.sorted.bam G_N2.sorted.bam H_N2.sorted.bam I_N2.sorted.bam
#/opt/samtools/x64/samtools/samtools index merged_bam_${run_ID}.sorted.bam

#eventually step to get rid of:
#/lscr2/andersenlab/kml436/test_Aug26/run_1/merged_bam_1.excision.cluster.rpmk.refined.bp.refsup
#/lscr2/andersenlab/kml436/test_Aug26/run_1/merged_bam_1.absence.refined.bp.summary
#/lscr2/andersenlab/kml436/test_Aug26/run_1/merged_bam_1.excision.cluster.rpmk.sfcp
#lscr2/andersenlab/kml436/test_Aug26/run_1/merged_bam_1.excision.cluster.rpmk.refined.bp
#/lscr2/andersenlab/kml436/test_Aug26/run_1/merged_bam_1.excision.cluster.rpmk

#mv merged_bam_${run_ID}.sorted.bam copy_merged_bam_${run_ID}.sorted.bam
#mv merged_bam_${run_ID}.sorted.bam.bai copy_merged_bam_${run_ID}.sorted.bam.bai
#java -jar /opt/picard-tools/DownsampleSam.jar P=.5 I=copy_merged_bam_${run_ID}.sorted.bam O=testing.bam
#/opt/samtools/x64/samtools/samtools sort -o -@ 8 testing.bam out > merged_bam_${run_ID}.sorted.bam
#/opt/samtools/x64/samtools/samtools index merged_bam_${run_ID}.sorted.bam

#create 2bit file of reference
#cd ${faToTwoBit_dir}
#faToTwoBit ${dir}/new_genome_${run_ID}.fasta new_genome_${run_ID}.2bit
#mv new_genome_${run_ID}.2bit ${dir}
#cd $dir

#remove preivous temp output so that it doesn't interfere with future steps
rm ${dir}/merged_bam_${run_ID}.excision.cluster.rpmk.refined.bp.refsup
rm ${dir}/merged_bam_${run_ID}.absence.refined.bp.summary
rm ${dir}/merged_bam_${run_ID}.excision.cluster.rpmk.sfcp
rm ${dir}/merged_bam_${run_ID}.excision.cluster.rpmk.refined.bp
rm ${dir}/merged_bam_${run_ID}.excision.cluster.rpmk


#run TEMP
echo "running TEMP..."
bash $absence -i merged_bam_${run_ID}.sorted.bam -s $TEMP_scripts -r adjusted_pos_all_tes_${run_ID}.bed -t new_genome_${run_ID}.2bit -c 20

#run TELOCATE
echo "running TELOCATE..."
echo "Creating sam file..."
mkdir sam_file
samtools view merged_bam_${run_ID}.sorted.bam|sort --temporary-directory=${dir}/sam_file > merged_sam_${run_ID}.sorted.sam
mkdir run_${run_ID}_filter_results_telocate/
mkdir run_${run_ID}_filter_results_temp/
mkdir run_${run_ID}_non_filter_results/
mv merged_sam_${run_ID}.sorted.sam sam_file
mkdir TELOCATE
cd ${TTR}/TE-locate/
perl TE_locate.pl 2 ${dir}/sam_file/ ${dir}/telocate_adjusted_pos_all_tes_${run_ID}.bed ${dir}/new_genome_${run_ID}.fasta ${dir}/TELOCATE/TEL $minimal_Distance_to_count $minimal_supporting_reads $minimal_supporting_individuals
cd ${dir}

cat merged_bam_${run_ID}.absence.refined.bp.summary | sed '1d' | awk -F $'\t' '{print $1"\t"$2"\t"$3"\t"$4"_reference_temp_rp_"NR"\t"$7"\t+"}' > run_${run_ID}_temp_nonredundant.bed
mv run_${run_ID}_temp_nonredundant.bed run_${run_ID}_filter_results_temp/
echo "TELOCATE finished...altering output..."
sed -e '1,2d' "${dir}/TELOCATE/TEL_${minimal_Distance_to_count}_reads${minimal_supporting_reads}_acc${minimal_supporting_individuals}.info" > "${dir}/run_${run_ID}_tmpfile"
awk -F'[\t/]' '{printf $1"\t"; if($17=="old") printf $2-1"\t"$2"\t"$8"\t"$5"_reference_telocate_rp_\t0"; else if($17=="new") printf $2-1"\t"$2"\t"$8"\t"$5"_non-reference_telocate_rp_\t0"; if($14=="parallel") print "\t+"; else if($14 =="uncertain") print "\t."; else print "\t-";}' "${dir}/run_${run_ID}_tmpfile" > "${dir}/run_${run_ID}_tmpfile_telocate_presort.txt"
bedtools sort -i "${dir}/run_${run_ID}_tmpfile_telocate_presort.txt" > "${dir}/run_${run_ID}_tmpfile_telocate_sorted_redundant.txt"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7}' "${dir}/run_${run_ID}_tmpfile_telocate_sorted_redundant.txt" > "${dir}/run_${run_ID}_tmpfile_telocate_counted_redundant.txt"
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' "${dir}/run_${run_ID}_tmpfile_telocate_counted_redundant.txt" >> "${dir}/run_${run_ID}_tmpfile_telocate_redundant.bed"
cat run_${run_ID}_tmpfile_telocate_redundant.bed | awk '$4~/_reference/ {print $0}' > tmp && mv tmp run_${run_ID}_tmpfile_telocate_redundant.bed
mv ${dir}/run_${run_ID}_tmpfile_telocate_redundant.bed run_${run_ID}_non_filter_results/
cd run_${run_ID}_non_filter_results/

#rename TELOCATE NON FILTER results (from WB element names to consensus names)....won't actualy end up using the "unfiltered" results...same as "filtered"
mv run_${run_ID}_tmpfile_telocate_redundant.bed run_${run_ID}_telocate_nonredundant.bed
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames run_${run_ID}_telocate_nonredundant.bed  TMP
mv run_${run_ID}_telocate_nonredundant.bed  ${run_ID}_nonredundant_org_copy.bed
mv TMP run_${run_ID}_telocate_nonredundant.bed
cat run_${run_ID}_telocate_nonredundant.bed| awk '$4 ~/_reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${run_ID}_temp && mv ${run_ID}_temp run_${run_ID}_telocate_nonredundant.bed

cd ..


echo "Filtering TELOCATE output......"
sort -k1,3 -k4rn "${dir}/run_${run_ID}_tmpfile_telocate_counted_redundant.txt"| sort -u -k1,3 | awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' > "${dir}/run_${run_ID}_tmpfileNR"
bedtools sort -i "${dir}/run_${run_ID}_tmpfileNR" >> "${dir}/run_${run_ID}_telocate_nonredundant.bed"
cat ${dir}/run_${run_ID}_telocate_nonredundant.bed | awk '$4~/_reference/ {print $0}' > tmp && mv tmp ${dir}/run_${run_ID}_telocate_nonredundant.bed
mv ${dir}/run_${run_ID}_telocate_nonredundant.bed run_${run_ID}_filter_results_telocate/


#rename TELOCATE results (from WB element names to consensus names)
cd run_${run_ID}_filter_results_telocate/
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames run_${run_ID}_telocate_nonredundant.bed  TMP
mv run_${run_ID}_telocate_nonredundant.bed  ${run_ID}_nonredundant_org_copy.bed
mv TMP run_${run_ID}_telocate_nonredundant.bed

#find closest matches to reference TEs within 1000 base pairs
closestBed -a ../adjusted_pos_all_tes_${run_ID}.bed -b run_${run_ID}_telocate_nonredundant.bed -d -t all | awk '$13<=1000 {print $0}' > ${run_ID}_closest.txt
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/remove_telocate_redundancies.py ${run_ID}_closest.txt
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/process_telocate_output.py closest_prepped.txt
cat processed_calls.txt| sort -k1,1 -k2,2n > run_${run_ID}_telocate_nonredundant.bed


cat run_${run_ID}_telocate_nonredundant.bed| awk '$4 ~/_reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}' > ${run_ID}_temp && mv ${run_ID}_temp run_${run_ID}_telocate_nonredundant.bed
cd ..

#rename TEMP results (from WB element names to consensus names)
cd run_${run_ID}_filter_results_temp/
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames run_${run_ID}_temp_nonredundant.bed  TMP
mv run_${run_ID}_temp_nonredundant.bed  ${run_ID}_nonredundant_org_copy.bed
mv TMP run_${run_ID}_temp_nonredundant.bed
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/process_double_deletion.py run_${run_ID}_temp_nonredundant.bed ../adjusted_pos_all_tes_${run_ID}.bed
mv tmp_double_deletion.txt run_${run_ID}_temp_nonredundant.bed
cat run_${run_ID}_temp_nonredundant.bed | awk '$4 ~/_reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}' > ${run_ID}_temp && mv ${run_ID}_temp run_${run_ID}_temp_nonredundant.bed

cd ..


# get the adjusted positions of the original reference TEs
cat adjusted_pos_all_tes_${run_ID}.bed | awk '$5=="R" && $4 !="TC8" {print $0}' > tmp && mv tmp adjusted_pos_REF_tes_${run_ID}.bed

# get the adjusted positions of the simulated reference TEs
cat adjusted_pos_all_tes_${run_ID}.bed | awk '$5!="R" && $4 !="TC8" {print $0}' > tmp && mv tmp adjusted_pos_SIM_tes_${run_ID}.bed
cp adjusted_pos_SIM_tes_${run_ID}.bed run_${run_ID}_filter_results_temp/
cp adjusted_pos_REF_tes_${run_ID}.bed run_${run_ID}_filter_results_telocate/

#FEED JUST SIMULATED to TEMP
cd run_${run_ID}_filter_results_temp/
python /lscr2/andersenlab/kml436/git_repos2/bamsurgeon/bed_compare9_redo_EDIT_RECALCULATE_absence.py 20 ${dir}/run_${run_ID}_filter_results_temp/ $TE_lengths adjusted_pos_SIM_tes_${run_ID}.bed run_${run_ID} $consensus_renamed >& run_${run_ID}_bedcompare.log
cd ..


####python /lscr2/andersenlab/kml436/git_repos2/bamsurgeon/bed_compare9_redo_EDIT_RECALCULATE_absence.py 20 /lscr2/andersenlab/kml436/RSV_TEST4/run_7/run_7_filter_results_temp/ /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt adjusted_pos_SIM_tes_7.bed run_7 /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/consensus.fasta >& run_7_bedcompare.log

#FEED JUST REFERENCE to TELOCATE
cp adjusted_pos_all_tes_${run_ID}.bed run_${run_ID}_filter_results_telocate/
cd run_${run_ID}_filter_results_telocate/
python /lscr2/andersenlab/kml436/git_repos2/bamsurgeon/bed_compare9_redo_EDIT_RECALCULATE_absence.py 20 ${dir}/run_${run_ID}_filter_results_telocate/ $TE_lengths adjusted_pos_REF_tes_${run_ID}.bed run_${run_ID} $consensus_renamed >& run_${run_ID}_bedcompare.log
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
cat run_${run_ID}_filter_results_temp/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed run_${run_ID}_filter_results_telocate/BEDCOMPARE_files_RECALCULATED/tmp_filter_file.bed > BEDCOMPARE_SUMMARY_run_${run_ID}_N2.txt
mv BEDCOMPARE_SUMMARY_run_${run_ID}_N2.txt final_results_RECALCULATED
cd final_results_RECALCULATED
cat FAMILY_TFPNs_F_temp FAMILY_TFPNs_F_telocate > FAMILY_TFPNs_F
cd ..


#for i in {1..1}
#do
#	echo $i
#	#sbatch
#	python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/WB_SIM_round5_TransposonBS_N2_SB_TTR.sh $i

#done
#bwa mem -t 24 new_genome_${run_ID}.fasta $fasta_1 $fasta_2 > test_${run_ID}.sam
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp test_${run_ID}.sam > test_${run_ID}.sorted.bam
#/opt/samtools/x64/samtools/samtools index test_${run_ID}.sorted.bam
#| /opt/samtools/x64/samtools/samtools view -@ 24 -bhu > test_1.unsorted.bam
###SEAPAERTE TEH ABOVE INTO 2 STEPS!!!!!
#/opt/samtools/x64/samtools/samtools sort -@ 24 -O bam -T tmp test_1.unsorted.bam > test_1.sorted.bam
#/opt/samtools/x64/samtools/samtools index test_1.sorted.bam
#/opt/samtools/x64/samtools/samtools sort -o -@ 8 test_1.unsorted.bam  out > test_1.sorted.bam
#/opt/samtools/x64/samtools/samtools index test_1.sorted.bam
