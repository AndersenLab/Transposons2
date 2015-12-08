#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=40000




##create error log file for everything?
### CHANGING JU11398 to JT11398

#######################
# FILE LOCATIONS
#######################
python_version=/exports/people/andersenlab/kml436/python2/bin/python
bam_surgeon=/lscr2/andersenlab/kml436/git_repos2/bamsurgeon
reference=/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa
twobit=/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.2bit
#twobit=/lscr2/andersenlab/kml436/11_c_elegans_reference.2bit #####
#2bit=/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.2bit
#2bit=/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.2bit 
TE_consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
processors=18
TTR=/lscr2/andersenlab/kml436/git_repos2/mcclintock/
HL_gff=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names_alias.bed
original_ref_pos=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names.gff
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt
family_renames=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt

minimal_Distance_to_count=1000
minimal_supporting_reads=3
minimal_supporting_individuals=1

#######################
# PREPARE BAM NAMING
#######################
bam_name=${1}
echo "bam name is ${bam_name}"
dir_to_bam=/lscr2/andersenlab/dec211/RUN/v2_snpset/bam
bam=${dir_to_bam}/${bam_name}.bam
echo "full path is ${bam}"
mkdir ${bam_name}
cd ${bam_name}
dir=`pwd`

cp -s $bam ${bam_name}.sorted.bam 
cp -s $bam.bai ${bam_name}.sorted.bam.bai

#######################
# RUN TEMP
#######################
mkdir raw_results_temp_insertion
cd raw_results_temp_insertion
echo "Running TEMP insertion caller..."
bash ${TTR}/TEMP/scripts/TEMP_Insertion.sh -i ../${bam_name}.sorted.bam -s ${TTR}/TEMP/scripts/ -x 30 -r $TE_consensus -t $HL_gff -m 1 -c 20 &> ../${bam_name}_TEMP_insertion_log.txt
cd ..

mkdir raw_results_temp_absence
cd raw_results_temp_absence
echo "Running TEMP absence caller..."
bash ${TTR}/TEMP/scripts/TEMP_Absence.sh -i ../${bam_name}.sorted.bam -s ${TTR}/TEMP/scripts/ -r $original_ref_pos -t $twobit -c 20 &> ../${bam_name}_TEMP_absence_log.txt
cd ..
#######################
# RUN TELOCATE
#######################
mkdir raw_results_telocate
cd raw_results_telocate
echo "Converting Bam to Sam..."
mkdir sam_file
cd sam_file
samtools view ../../${bam_name}.sorted.bam |sort --temporary-directory=${dir}/raw_results_telocate/sam_file > ${bam_name}.sorted.sam
cd ..
echo "Running TELOCATE absence caller..."
mkdir TELOCATE
cd ${TTR}/TE-locate/
perl TE_locate.pl 2 ${dir}/raw_results_telocate/sam_file/ $HL_gff $reference ${dir}/raw_results_telocate/TELOCATE/TEL $minimal_Distance_to_count $minimal_supporting_reads $minimal_supporting_individuals &> ${dir}/${bam_name}_TELOCATE_log.txt
cd ${dir}


#remove sam files and copied-over bam files to save space
###CHECK THE BELOW
rm raw_results_temp_insertion/${bam_name}.sorted.bam
rm raw_results_temp_insertion/${bam_name}.sorted.bam.bai
rm raw_results_temp_absence/${bam_name}.sorted.bam
rm raw_results_temp_absence/${bam_name}.sorted.bam.bai
########UNDO THE BELOW LINE ALTER!!!!!!!#########
rm -rf raw_results_telocate/sam_file

echo "Done"

#ADD IN TELOCATE
#ADD IN TEMP EXCISION
#mkdir final_results
