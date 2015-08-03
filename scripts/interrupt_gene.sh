#!/bin/bash
# this script checks which genes have transposon inserted or removed from them
# outputs each method separately and output one concatenated file as well
# USE: interrupt_gene.sh <sample_list>

file=${1}
dir=/lscr2/andersenlab/kml436/git_repos2/Transposons/data
filedir=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/
gene_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/gene.gff
promoter_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/promoter.gff
intron_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/intron.gff
exon_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/exon.gff
piRNA_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/piRNA.gff
process_script=/lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/process_gene_interrupt.py

#pull only Wormbase marked positions!!!
cat  $gene_positions |awk '$2=="WormBase" {print $0}' > ${filedir}WB_gene_positions.gff
cat  $promoter_positions > ${filedir}WB_promoter_positions.gff
cat  $intron_positions|awk '$2=="WormBase" {print $0}' > ${filedir}WB_intron_positions.gff
cat  $exon_positions |awk '$2=="WormBase" {print $0}' > ${filedir}WB_exon_positions.gff
cat  $piRNA_positions |awk '$2=="WormBase" {print $0}' > ${filedir}WB_piRNA_positions.gff

#remove existing results  files....change this later
rm all_gene_insertions_and_excisions.txt
rm all_promoter_insertions_and_excisions.txt
rm all_introns_insertions_and_excisions.txt
rm all_exons_insertions_and_excisions.txt
rm all_piRNAs_insertions_and_excisions.txt


sed -e '/^#/d' -e '/^$/d' $file > tmp && mv tmp $file


# simplify with loop for everything but gene depending on wanted window
while read -r sample
do
	#GENE FILE
	# use bedtools window with -w to extend search by certain number of base pairs on either side
	bedtools window -a ${filedir}WB_gene_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_insertion_nonredundant.bed -l 1000 -r 25 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tinsertion\t"$0}' > insertions_into_genes.txt
	bedtools window -a ${filedir}WB_gene_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_absence_nonredundant.bed -w 25 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tabsence\t"$0}' > absences_from_genes.txt
	bedtools window -a ${filedir}WB_gene_positions.gff -b ${dir}/${sample}/final_results/${sample}_telocate_nonredundant.bed -w 25 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\treference\t"$0}' > references_in_genes.txt

	#PROMOTER FILE
	bedtools window -a ${filedir}WB_promoter_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_insertion_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tinsertion\t"$0}' > insertions_into_promoters.txt
	bedtools window -a ${filedir}WB_promoter_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_absence_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tabsence\t"$0}' > absences_from_promoters.txt
	bedtools window -a ${filedir}WB_promoter_positions.gff -b ${dir}/${sample}/final_results/${sample}_telocate_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\treference\t"$0}' > references_in_promoters.txt

	#INTRON FILE
	bedtools window -a ${filedir}WB_intron_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_insertion_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tinsertion\t"$0}' > insertions_into_introns.txt
	bedtools window -a ${filedir}WB_intron_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_absence_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tabsence\t"$0}' > absences_from_introns.txt
	bedtools window -a ${filedir}WB_intron_positions.gff -b ${dir}/${sample}/final_results/${sample}_telocate_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\treference\t"$0}' > references_in_introns.txt


	#EXON FILE
	bedtools window -a ${filedir}WB_exon_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_insertion_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tinsertion\t"$0}' > insertions_into_exons.txt
	bedtools window -a ${filedir}WB_exon_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_absence_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tabsence\t"$0}' > absences_from_exons.txt
	bedtools window -a ${filedir}WB_exon_positions.gff -b ${dir}/${sample}/final_results/${sample}_telocate_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\treference\t"$0}' > references_in_exons.txt

	#piRNAFILE
	bedtools window -a ${filedir}WB_piRNA_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_insertion_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tinsertion\t"$0}' > insertions_into_piRNAs.txt
	bedtools window -a ${filedir}WB_piRNA_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_absence_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tabsence\t"$0}' > absences_from_piRNAs.txt
	bedtools window -a ${filedir}WB_piRNA_positions.gff -b ${dir}/${sample}/final_results/${sample}_telocate_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\treference\t"$0}' > references_in_piRNAs.txt

	#put together
	cat insertions_into_genes.txt absences_from_genes.txt references_in_genes.txt | cut -f1-3,5,11,13,15 >> all_gene_insertions_and_excisions.txt
	cat insertions_into_promoters.txt absences_from_promoters.txt references_in_promoters.txt | cut -f1-3,5,11,13,15 >> all_promoter_insertions_and_excisions.txt
	cat insertions_into_introns.txt absences_from_introns.txt references_in_introns.txt | cut -f1-3,5,11,13,15 >> all_introns_insertions_and_excisions.txt
	cat insertions_into_exons.txt absences_from_exons.txt references_in_exons.txt | cut -f1-3,5,11,13,15 >> all_exons_insertions_and_excisions.txt
	cat insertions_into_piRNAs.txt absences_from_piRNAs.txt references_in_piRNAs.txt | cut -f1-3,5,11,13,15 >> all_piRNAs_insertions_and_excisions.txt

done < "${file}"

python $process_script all_gene_insertions_and_excisions.txt
python $process_script all_promoter_insertions_and_excisions.txt
python $process_script all_introns_insertions_and_excisions.txt
python $process_script all_exons_insertions_and_excisions.txt
python $process_script all_piRNAs_insertions_and_excisions.txt



cat final_all_gene_insertions_and_excisions.txt\
 final_all_promoter_insertions_and_excisions.txt\
  final_all_introns_insertions_and_excisions.txt\
   final_all_exons_insertions_and_excisions.txt\
   final_all_piRNAs_insertions_and_excisions.txt > final_gene_table.txt