#!/bin/bash
# this script runs the ProcessTransposonCallerOutput script on each sample specified in the input_file and outputs the summary output to the specified results_file

input_file=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/test_list.txt
ProcessTransposonCallerOutput=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/alt_ProcessTransposonCallerOutput.sh
results_file=/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/FINAL_RESULTS.txt



### avoid blank lines and comments
sed -e '/^#/d' -e '/^$/d' $input_file > tmp && mv tmp $input_file
while read line; do 
  echo $line
  echo "Running ProcessTransposonCaller on $line"
  sbatch $ProcessTransposonCallerOutput $line $results_file
done <$input_file
