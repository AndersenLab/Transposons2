#!/bin/bash
# this script calculates the averages of multiple simulation runs (such as average TPR and FDR)
# and genreates graphs/plots of the results with R
# NOTE: comment out "calculate_mean("retroseq")" and "calculate_mean("telocate")" in family_TFPN_average4.py and TFPN_average4.py as needed
# NOTE: family_average script doesn't always require using the renamed consensus seqs
# USE:analyze_simulations_ROUND18.sh in appropriate folder

#recalculate_edit=RECALCULATE_EDIT.sh
#recalculate=RECALCULATE.sh
average=run_average_TFPNs_ROUND18.sh
family_average=run_family_TFPNs_ROUND18.sh
graph_TFPN=graph_TFPR_distances_ROUND18.R
graph_family_TFPN=graph_TFPN_families.R
script_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons/scripts


consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/round2_consensus_set2.fasta
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/AB-PR/consensus_wTC8.fasta 
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/AB-PR/fake_lengths.txt

#bash ${script_dir}/${recalculate} $length
#merge the distance calculations
mkdir averages
#calculate averages
bash ${script_dir}/${average} $consensus # consensus seq file doesn't matter here
#only sometimes run the below line
bash ${script_dir}/${family_average} $consensus_renamed
#graph the TPR and FDR
cd averages
Rscript ${script_dir}/${graph_TFPN}
cd ..
#plot the TPR and FDR per family
#only sometimes run the below line
#cd families
#Rscript ${script_dir}/${graph_family_TFPN}
#cd ..

