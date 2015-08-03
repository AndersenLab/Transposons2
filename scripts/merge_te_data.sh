#!/bin/bash
# this script merges the nonredundant te files into one file and adds a columns to distinguish the caller method associated with the data
# USE: merge_te_data.sh in appropriate directory

find ./*/final_results/ -name "*_temp_insertion_nonredundant.bed" -exec awk '{print $0"\tnew"}' {} > intermediate.txt \;
find ./*/final_results/ -name "*_temp_absence_nonredundant.bed" -exec awk '{print $0"\tabsent"}' {} >> intermediate.txt \;
find ./*/final_results/ -name "*_telocate_nonredundant.bed" -exec awk '{print $0"\treference"}' {} >> intermediate.txt \;
mv intermediate.txt all_nonredundant.txt