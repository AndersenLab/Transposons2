#!/bin/bash
# this script moves the neccesary files from each run into a "families" folder and runs the family_TFPN_average4.py script on these files
# USE: run_family_TFPNS.sh <consensus_file>
# called on by analyze_simulations.sh

consensus=${1}
#consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/round3_consensus_fasta.fasta
#consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/round2_consensus_fasta.fasta
mkdir families
for i in {1..8}
do
	cp run_${i}_N2/final_results_RECALCULATED/FAMILY_TFPNs_F families/FAMILY_TFPN_ALL_run_${i}_N2

done
cd families
dir=`pwd`
echo $dir
python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/family_TFPN_average4.py $dir $consensus
