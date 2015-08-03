#!/bin/bash
# this script runs the transposon simulation and detection scripts with sbatch the number of times specified (edit script)
# usage: run_simulations.sh
for i in {1..10}
do
	echo $i
	#sbatch /exports/people/andersenlab/kml436/scripts/TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/round2_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/round3_TransposonBS_N2_SB_TTR.sh $i
	sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round2_TransposonBS_N2_SB_TTR.sh $i

done
