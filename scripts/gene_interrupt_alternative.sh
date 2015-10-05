#!/bin/bash
# this script runs the kin matrix script and checks if the unique TEs are in genic or intergenic regions
# USE: gene_interrupt_alternative.sh


filedir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/
script_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/
kin_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/kintest/TE_matrix/


cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/kintest
#in the kin_mean.py script set the distance equal to 10!!!!
python ${script_dir}kin_mean.py ${kin_dir}/../CtCp_all_nonredundant.txt

#produce bedfile for insertion, absence, and reference calls
cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/kintest/TE_matrix
python ${script_dir}split_bedt.py
cat bedfile_insertions.bed bedfile_absences.bed bedfile_references.bed > all.bed
cd ..
cp TE_matrix/all.bed .

# pull only Wormbase protein_coding and psedogenes for gene gff file:
cat ${filedir}gene.gff | awk '$2=="WormBase" && $9 ~/biotype=(protein_coding|pseudogene|transposon_pseudogene)/  {print $0}'> ${filedir}WB_gene_positions.gff

# created clipped gene file contained the original base postions clipped by 10 on either side
cd $filedir
python ${script_dir}modify_gff.py

cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/kintest
# run bedtools with window set to 0 for absences and references, 0 or 10 for insertions, one for full one for clipped
#absences
bedtools window -a ${filedir}WB_gene_positions.gff -b ${kin_dir}bedfile_absences.bed -w 0 > absences_window.txt
#references
bedtools window -a ${filedir}WB_gene_positions.gff -b ${kin_dir}bedfile_references.bed -w 0 > references_window.txt
#insertions
bedtools window -a ${filedir}clipped_WB_gene_positions.gff -b ${kin_dir}bedfile_insertions.bed -w 0 > insertions_window.txt
bedtools window -a ${filedir}clipped_WB_gene_positions.gff -b ${kin_dir}bedfile_insertions.bed -w 10 > insertions_full_window.txt

cat absences_window.txt references_window.txt insertions_window.txt > all_window.txt

#classify whether a TE is genic or intergenic, and if geneic if internal or border(wihtin 10 bps of the end)
python ${script_dir}separate_gene_assignments.py

