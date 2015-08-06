#!/bin/bash
# this script merges the nonredundant te files into one file and adds a columns to distinguish the caller method associated with the data
# USE: merge_te_data.sh in appropriate directory (usually the data directory)

find ./*/final_results/ -name "*_temp_insertion_nonredundant.bed" -exec awk -v path={}  '{match(path,"\.\/(.*)\/final_results\/.*",a)}; {print $0"\tnew\t"a[1]}' {} > intermediate.txt \;
find ./*/final_results/ -name "*_temp_absence_nonredundant.bed" -exec awk -v path={}  '{match(path,"\.\/(.*)\/final_results\/.*",a)}; {print $0"\tabsent\t"a[1]}' {} >> intermediate.txt \;
find ./*/final_results/ -name "*_telocate_nonredundant.bed" -exec awk -v path={}  '{match(path,"\.\/(.*)\/final_results\/.*",a)};{print $0"\treference\t"a[1]}' {} >> intermediate.txt \;
mv intermediate.txt all_nonredundant.txt


#alternatives
#find ./*/final_results/ -name "*_temp_insertion_nonredundant.bed" -exec awk '{print $0"\tnew"}' {} > intermediate.txt \;
#find ./*/final_results/ -name "*_temp_absence_nonredundant.bed" -exec awk '{print $0"\tabsent"}' {} >> intermediate.txt \;
#find ./*/final_results/ -name "*_telocate_nonredundant.bed" -exec awk '{print $0"\treference"}' {} >> intermediate.txt \;
#find ./*/final_results/ -name "*_temp_insertion_nonredundant.bed" -exec echo {} \;  -exec awk -v path={}  '{match(path,"\.\/(.*)\/final_results\/.*",a)}END{print a[1]}' {} \;
#find ./*/final_results/ -name "*_temp_insertion_nonredundant.bed" -exec echo {} \;  -exec awk -v path={} '{print $0"\tnew\t"path}'  {} \;


