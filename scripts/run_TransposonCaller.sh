#!/bin/bash
# this script runs the TransposonCaller script on each sample specified in the input_file
# USE: run_TransposonCaller.sh in the data dir

input_file=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/master_sample_list.txt
TransposonCaller=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/TransposonCaller.sh



### avoid blank lines and comments
sed -e '/^#/d' -e '/^$/d' $input_file > tmp && mv tmp $input_file
while read line; do 
  echo $line
  echo "Running TransposonCaller on $line"
  sbatch $TransposonCaller $line
done <$input_file



