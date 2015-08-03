#!/bin/bash
# this script moves the neccesary files from each run into a "averages" folder and runs the TFPN_average4.py script on these files
# USE: run_family_TFPNS.sh <consensus_file>
# called on by analyze_simulations_ROUND18.sh
# consensus file doesn't matter here
# NOTE: specific to ROUND18 data output

consensus=${1}
#consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/round2_consensus_fasta.fasta
#consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/round3_consensus_fasta.fasta
#consensus=WB_all_tes_plus_strand.fasta
mkdir averages
for i in {1..1}
do
	cp run_${i}_N2/final_results_RECALCULATED/BEDCOMPARE_SUMMARY_run_${i}_N2.txt averages/

done
cd averages
dir=`pwd`
echo $dir
python /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/TFPN_average4.py $dir 
