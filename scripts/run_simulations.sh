#!/bin/bash
# this script runs the transposon simulation and detection scripts with sbatch the number of times specified (edit script)
# USE: run_simulations.sh
for i in {1..8}
do
	echo $i
	#sbatch /exports/people/andersenlab/kml436/scripts/TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/round2_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/round3_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round4_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round5_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round6_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round7_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round8_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round9_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round10_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round11_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round12_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round13_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round14_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round15_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round16_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round17_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/run_RSVSim.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/RSV_ReadSupport_Loop.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/WB_SIM_round18_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/ROUND18_ReadSupport_Loop.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/WB_SIM_round19_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/WB_SIM_round20_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/RSV_ReadSupport_Loop.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/RSV_VariantSupport_Loop_round20.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/WB_SIM_round21_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/WB_SIM_round22_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/WB_SIM_round23_TransposonBS_N2_SB_TTR.sh $i
	#sbatch /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/RSV_ReadSupport_Loop_round20.sh $i

done
