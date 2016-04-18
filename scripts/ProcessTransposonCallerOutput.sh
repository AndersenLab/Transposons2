#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --nodelist=node6
#SBATCH --mem=2000

family_renames=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 
#consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/consensus.fasta
ref_tes=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_family_names.gff

minimal_Distance_to_count=1000
minimal_supporting_reads=3
minimal_supporting_individuals=1

bam_name=${1}
results_file=${2}

cd $bam_name
mkdir final_results
cd raw_results_temp_insertion


##############################################################################################
# 						PROCESSING TEMP INSERTION OUTPUT
##############################################################################################

echo "Processing TEMP output..."
rm ${bam_name}_temp_nonredundant.bed
rm ${bam_name}_temp_raw.bed
rm ${bam_name}_temp_redundant.bed

cp ${bam_name}.insertion.refined.bp.summary current_${bam_name}.insertion.refined.bp.summary
sed '1d' "current_${bam_name}.insertion.refined.bp.summary" | awk '{if ($4 == "sense" || $4 == "antisense"); else print $0}' | awk '{ printf $1"\t"$9"\t"$11"\t"$7"\t"$4"_non-reference_"sample"_temp_\t0\t"; if ( $5 == "sense" ) printf "+"; else printf "-"; print "\t"$6"\t"$10"\t"$12"\t"$8}' > "${bam_name}_temp_presort_raw.txt"
bedtools sort -i "${bam_name}_temp_presort_raw.txt" > "${bam_name}_temp_sorted_raw.txt"
awk '{ printf $1"\t"$2"\t"$3"\t"$4"\t"$5; if ($9 > 0 && $10 > 0) printf "sr_"; else printf "rp_"; print NR"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' "${bam_name}_temp_sorted_raw.txt" > "${bam_name}_temp_sorted_precut_raw.txt"
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7"\t"$11}' "${bam_name}_temp_sorted_precut_raw.txt" >> "${bam_name}_temp_raw.bed"

echo "Filtering TEMP output......"
# Filter for insertions with support for both ends
awk '{if ($8 == "1p1") print $0}' "${bam_name}_temp_sorted_precut_raw.txt" > "${bam_name}_temp_sorted_redundant.txt"
cut -f1-3,5-7 "${bam_name}_temp_sorted_redundant.txt" >> "${bam_name}_temp_redundant.bed"
# Filter out redundant insertions
sort -k1,1 -k2,2n -k4rn "${bam_name}_temp_sorted_redundant.txt" | sort -u -k1,2| awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7"\t"$11}' > "${bam_name}_tmp"
bedtools sort -i "${bam_name}_tmp">> "${bam_name}_temp_nonredundant.bed"


# RENAME TRANSPOSON CALLS
#rename the TEs in the output tp match final consensus names
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${bam_name}_temp_nonredundant.bed TMP
mv ${bam_name}_temp_nonredundant.bed ${bam_name}_nonredundant_org_copy.bed
mv TMP ${bam_name}_temp_nonredundant.bed
#rename chr X, ensure all lines contain "non-reference" (should already be the case for TEMP)
cat ${bam_name}_temp_nonredundant.bed| sed s'/chrX/X/g'| sed '/^track/d' | awk '$4 ~/non-reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${bam_name}_temp && mv ${bam_name}_temp ${bam_name}_temp_nonredundant.bed  
cat ${bam_name}_temp_nonredundant.bed | awk '$7>.25 && $5>=8 {print $0}' > tmp && mv tmp ${bam_name}_temp_nonredundant.bed #rec_edit

# COLLAPSE TRANSPOSONS
#run collapse transposons script
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/collapse_transposons.py $TE_lengths ${bam_name}_temp_nonredundant.bed
cp new_CT_${bam_name}_temp_nonredundant.bed ../final_results/${bam_name}_temp_insertion_nonredundant.bed

# COUNT TRANSPOSONS
#run te_totals.py script
end_tes=`grep -c "specifies an unknown reference name. Continue anyway." ../${bam_name}_TEMP_insertion_log.txt`
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/te_totals.py new_CT_${bam_name}_temp_nonredundant.bed $consensus_renamed ${bam_name} $results_file new $end_tes
cd ..


##############################################################################################
# 						PROCESSING TEMP ABSENCE OUTPUT
##############################################################################################
cd raw_results_temp_absence
rm ${bam_name}_temp_nonredundant.bed
cat ${bam_name}.absence.refined.bp.summary | sed '1d' | awk -F $'\t' '{print $1"\t"$2"\t"$3"\t"$4"_reference_temp_rp_"NR"\t"$7"\t+"}' > ${bam_name}_temp_nonredundant.bed 

# RENAME TRANSPOSON CALLS
#rename the TEs in the output tp match final consensus names
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${bam_name}_temp_nonredundant.bed TMP
mv ${bam_name}_temp_nonredundant.bed ${bam_name}_nonredundant_org_copy.bed
mv TMP ${bam_name}_temp_nonredundant.bed
#rename chr X, ensure all lines contain "reference" 
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/process_double_deletion.py ${bam_name}_temp_nonredundant.bed $ref_tes
mv tmp_double_deletion.txt ${bam_name}_temp_nonredundant.bed
cat ${bam_name}_temp_nonredundant.bed| sed s'/chrX/X/g'| sed '/^track/d' | awk '$4 ~/_reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${bam_name}_temp && mv ${bam_name}_temp ${bam_name}_temp_nonredundant.bed 
cat ${bam_name}_temp_nonredundant.bed | sort -k1,1 -k2,2n > ${bam_name}_temp && mv ${bam_name}_temp ${bam_name}_temp_nonredundant.bed 

# Filter out redundant calls (which can occur if a TE was counted both indiviually and in a 'double deletion' call, read  support is in column 5 )
sort -k1,1 -k2,2n -k5rn "${bam_name}_temp_nonredundant.bed" | sort -u -k1,2 > "${bam_name}_tmp"
bedtools sort -i "${bam_name}_tmp"> tmp && mv tmp ${bam_name}_temp_nonredundant.bed


cat ${bam_name}_temp_nonredundant.bed | awk '$5>=3 {print $0}' > tmp && mv tmp ${bam_name}_temp_nonredundant.bed #rec_edit
cp ${bam_name}_temp_nonredundant.bed ../final_results/${bam_name}_temp_absence_nonredundant.bed


# COUNT TRANSPOSONS
#run te_totals.py script
end_tes=`grep -c "specifies an unknown reference name. Continue anyway." ../${bam_name}_TEMP_absence_log.txt`
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/te_totals.py ${bam_name}_temp_nonredundant.bed $consensus_renamed ${bam_name} $results_file absent $end_tes
cd ..

##############################################################################################
# 						PROCESSING TELOCATE OUTPUT
##############################################################################################

cd raw_results_telocate
rm ${bam_name}_telocate_nonredundant.bed
rm ${bam_name}_tmpfile_telocate_redundant.bed
echo "Filtering TELOCATE output..."
sed -e '1,2d' "./TELOCATE/TEL_${minimal_Distance_to_count}_reads${minimal_supporting_reads}_acc${minimal_supporting_individuals}.info" > "${bam_name}_tmpfile"
awk -F'[\t/]' '{printf $1"\t"; if($17=="old") printf $2-1"\t"$2"\t"$8"\t"$5"_reference_telocate_rp_\t0"; else if($17=="new") printf $2-1"\t"$2"\t"$8"\t"$5"_non-reference_telocate_rp_\t0"; if($14=="parallel") print "\t+"; else if($14 =="uncertain") print "\t."; else print "\t-";}' "${bam_name}_tmpfile" > "${bam_name}_tmpfile_telocate_presort.txt"
bedtools sort -i "${bam_name}_tmpfile_telocate_presort.txt" > "${bam_name}_tmpfile_telocate_sorted_redundant.txt"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7}' "${bam_name}_tmpfile_telocate_sorted_redundant.txt" > "${bam_name}_tmpfile_telocate_counted_redundant.txt"
awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}' "${bam_name}_tmpfile_telocate_counted_redundant.txt" >> "${bam_name}_tmpfile_telocate_redundant.bed"
cat ${bam_name}_tmpfile_telocate_redundant.bed | awk '$4~/_reference/ {print $0}' > tmp && mv tmp ${bam_name}_tmpfile_telocate_redundant.bed 
sort -k1,1 -k2,2n -k4rn  "${bam_name}_tmpfile_telocate_counted_redundant.txt"| sort -u -k1,2| awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$7}'> "${bam_name}_tmpfileNR"
bedtools sort -i "${bam_name}_tmpfileNR" >> "${bam_name}_telocate_nonredundant.bed"
cat ${bam_name}_telocate_nonredundant.bed | awk '$4~/_reference/ {print $0}' > tmp && mv tmp ${bam_name}_telocate_nonredundant.bed 

# RENAME TRANSPOSON CALLS
#rename the TEs in the output tp match final consensus names
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py $family_renames ${bam_name}_telocate_nonredundant.bed TMP
mv ${bam_name}_telocate_nonredundant.bed ${bam_name}_nonredundant_org_copy.bed 
mv TMP ${bam_name}_telocate_nonredundant.bed 
#rename chr X, ensure all lines contain "reference" 
closestBed -a $ref_tes -b ${bam_name}_telocate_nonredundant.bed  -d -t all | awk '$13<=1000 {print $0}' > ${bam_name}_closest.txt
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/process_telocate_output.py ${bam_name}_closest.txt
cat processed_calls.txt| sort -k1,1 -k2,2n > ${bam_name}_telocate_nonredundant.bed 


cat ${bam_name}_telocate_nonredundant.bed | sed s'/chrX/X/g'| sed '/^track/d' | awk '$4 ~/_reference/ && $4 !~ /TC8/ && $1 !~ /MtDNA/ {print $0}'> ${bam_name}_temp && mv ${bam_name}_temp ${bam_name}_telocate_nonredundant.bed 
cp ${bam_name}_telocate_nonredundant.bed ../final_results/${bam_name}_telocate_nonredundant.bed ###HERE ERROR


# COUNT TRANSPOSONS
#run te_totals.py script
end_tes=`grep -c "specifies an unknown reference name. Continue anyway." ../${bam_name}_TELOCATE_log.txt`
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/te_totals.py ${bam_name}_telocate_nonredundant.bed  $consensus_renamed ${bam_name} $results_file reference $end_tes
cd ..

echo "Done"
