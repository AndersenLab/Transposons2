#!/bin/bash
# old check script...ignore

run_ID=${1}
mkdir run_${run_ID}
cd run_${run_ID}
dir=`pwd`
echo $dir

distance=${2}


BEDCOMPARE_SUMMARY_run_${run_ID}_N2.txt 
# look at FAMILY_TFPNs_F later
cat BEDCOMPARE_SUMMARY_run_${run_ID}.txt | awk -v distance="$distance" '$3 >=distance {print $0}' > tmp_filter_file.bed