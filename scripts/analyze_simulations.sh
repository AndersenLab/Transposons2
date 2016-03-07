#!/bin/bash
# this script calculates the averages of multiple simulation runs (such as average TPR and FDR)
# and genreates graphs/plots of the results with R
# NOTE: comment out "calculate_mean("retroseq")" and "calculate_mean("telocate")" in family_TFPN_average4.py and TFPN_average4.py
# NOTE: family_average script doesn't always require using the renamed consensus seqs
# USE:analyze_simulations.sh in appropriate folder

recalculate_edit=RECALCULATE_EDIT_NEWCOLLAPSE.sh

#recalculate_edit=RECALCULATE_EDIT.sh
#recalculate=RECALCULATE.sh
merge=merge_SC_files.sh
average=run_average_TFPNs.sh
family_average=run_family_TFPNs.sh
graph_TFPN=graph_TFPR_distances.R
graph_family_TFPN=graph_TFPN_families.R
script_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts


#consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
#consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/consensus.fasta
#length=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 
length=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths_plusTC8.txt

#adjust positions so that only include the start and not the end of the detected TE 
bash ${script_dir}/${recalculate_edit} $consensus_renamed $length 100
#bash ${script_dir}/${recalculate} $length
#merge the distance calculations
bash ${script_dir}/${merge}
mkdir averages
cp filtered_SCs/SCs_ALL_filter averages
cp non_filtered_SCs/SCs_ALL_non_filter averages
#calculate averages
bash ${script_dir}/${average} $consensus # consensus seq file doesn't matter here
#only sometimes run the below line
#bash ${script_dir}/${family_average} $consensus_renamed
#graph the TPR and FDR
#cd averages
#Rscript ${script_dir}/${graph_TFPN}
#cd ..
#plot the TPR and FDR per family
#only sometimes run the below line
#cd families
#Rscript ${script_dir}/${graph_family_TFPN}
#cd ..

#bash /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/RECALCULATE_EDIT.sh /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths.txt
#bash /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/merge_SC_files.sh
#bash /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/run_average_TFPNs.sh /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
#bash /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/run_family_TFPNs.sh /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
#Rscript /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/graph_TFPN_families.R
