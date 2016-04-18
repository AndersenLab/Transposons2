#!/bin/bash
# this script analyzes GWAS results to investigate invlovement of genes of interest and piRNAs
# USE: PostMappings.sh
# NOTE: must first transfer "Peak_Table.txt" to final_results directory and "vc_PI.txt" to piRNA directory


scripts_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts

# check if QTL overlap regions of interest
echo "Investigating QTL overlaps..."
python ${scripts_dir}/QTL_control_overlap.py

# check if piRNA align to TE sequences
echo "Investigating piRNA alignments to TEs..."
cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA
mkdir all
cd all
python ${scripts_dir}/pi_v_TE.py

# check alignment of piRNA to respective TE families
echo "Investigating QTL specific piRNA alignments..."
cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA
python ${scripts_dir}/piBWA.py

# BLAST all piRNA to all TE seqs
echo "Investigating piRNA blasted to TE seqs..."
mkdir blast
cd blast
python ${scripts_dir}/pi_BLAST_all.py