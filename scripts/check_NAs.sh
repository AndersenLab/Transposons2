#!/usr/bin/bash
# this script checks which samples have either an absence or a reference call at a postion....if not, counted as "NA"
# USE: check_NAs.sh

file=${1}
data_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data
cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/results

touch non_Nas.txt 

while read -r sample
do
	cd $data_dir
	cd ${sample}/final_results
	ref_calls=`wc -l ${sample}_telocate_nonredundant.bed|cut -d' ' -f1`
	abs_calls=`wc -l ${sample}_temp_absence_nonredundant.bed|cut -d' ' -f1`
	total=$(echo "scale=4; $ref_calls + $abs_calls" | bc -l)
	
	overall="$sample $total"
	echo $overall >> /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/non_Nas.txt
done < "${file}"

cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/results
cat non_Nas.txt| sort -k2,2nr > tmp && mv tmp non_Nas.txt
