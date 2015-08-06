#!/bin/bash
# this script recount the transposon totals after running the script to resolve contradictory calls between the telocate reference caller and temp absence caller
# USE: run_recount.sh <sample_list file> <name of results file> <consensus file> <data directory>

input_file=${1} #/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/full_sample_list.txt
results_file=${2} #/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/FINAL_RESULTS.txt
consensus_renamed=${3} #/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 
data_dir=${4} #/lscr2/andersenlab/kml436/git_repos2/Transposons2/data

### avoid blank lines and comments
sed -e '/^#/d' -e '/^$/d' $input_file > tmp && mv tmp $input_file
while read line; do 
  echo "Counting Transposons in $line"
  #run te_totals.py script TEMP insertion
  end_tes=`grep -c "specifies an unknown reference name. Continue anyway." ${data_dir}/${line}/${line}_TEMP_insertion_log.txt`
  python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/te_totals.py ${data_dir}/${line}/final_results/${line}_temp_insertion_nonredundant.bed $consensus_renamed ${line} $results_file new $end_tes

  #run te_totals.py script TELOCATE reference
  end_tes=`grep -c "specifies an unknown reference name. Continue anyway." ${data_dir}/${line}/${line}_TELOCATE_log.txt`
  python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/te_totals.py ${data_dir}/${line}/final_results/${line}_telocate_nonredundant.bed  $consensus_renamed ${line} $results_file reference $end_tes

  #run te_totals.py script TEMP absence
  end_tes=`grep -c "specifies an unknown reference name. Continue anyway." ${data_dir}/${line}/${line}_TEMP_absence_log.txt`
  python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/te_totals.py ${data_dir}/${line}/final_results/${line}_temp_absence_nonredundant.bed $consensus_renamed ${line} $results_file absent $end_tes
done <$input_file

