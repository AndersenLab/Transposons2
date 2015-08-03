#!/bin/bash
# this script mereges the SC files from a simulation set (treats filtered and non-filtered results separately)
# USE: merge_SC_files.sh in appropriate directory

mkdir filtered_SCs
mkdir non_filtered_SCs
for i in {1..8}
do
	cp run_${i}_N2/run_${i}_N2_filter_results/intermediates_RECALCULATED/SC* filtered_SCs/
	cp run_${i}_N2/run_${i}_N2_non_filter_results/intermediates_RECALCULATED/SC* non_filtered_SCs/

done
cd filtered_SCs
cat SC* > SCs_ALL_filter
cd ..
cd non_filtered_SCs
cat SC* > SCs_ALL_non_filter
cd ..
