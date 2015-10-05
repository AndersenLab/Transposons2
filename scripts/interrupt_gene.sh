# this script checks which genes,introns, exons, five prime UTRs, and threeprime UTRs have transposon inserted or removed from them
# outputs each method separately and output one concatenated file as well
# USE: interrupt_gene.sh <sample_list>

file=${1}
dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data
filedir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/
gene_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/gene.gff
intron_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/intron.gff
exon_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/exon.gff
fiveUTR_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/five_prime_UTR.gff
threeUTR_positions=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/three_prime_UTR.gff
process_script=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/process_gene_interrupt.py
check_script=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/gene_check.py
other_gene=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/other_gene_determination.py

#pull only Wormbase marked positions!!!
cat  $gene_positions |awk '$2=="WormBase" {print $0}' > ${filedir}WB_gene_positions.gff
cat  $intron_positions|awk '$2=="WormBase" {print $0}' > ${filedir}WB_intron_positions.gff
cat  $exon_positions |awk '$2=="WormBase" {print $0}' > ${filedir}WB_exon_positions.gff
cat  $fiveUTR_positions |awk '$2=="WormBase" {print $0}' > ${filedir}WB_fiveUTR_positions.gff
cat  $threeUTR_positions |awk '$2=="WormBase" {print $0}' > ${filedir}WB_threeUTR_positions.gff


#remove existing results  files....change this later
rm all_gene_insertions_and_excisions.txt
rm all_introns_insertions_and_excisions.txt
rm all_exons_insertions_and_excisions.txt
rm all_fiveUTR_insertions_and_excisions.txt
rm all_threeUTR_insertions_and_excisions.txt


sed -e '/^#/d' -e '/^$/d' $file > tmp && mv tmp $file


# simplify with loop for everything but gene depending on wanted window
while read -r sample
do
	#GENE FILE
	# use bedtools window with -w to extend search by certain number of base pairs on either side
	bedtools window -a ${filedir}WB_gene_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_insertion_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tinsertion\t"$0}' > insertions_into_genes.txt
	bedtools window -a ${filedir}WB_gene_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_absence_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tabsence\t"$0}' > absences_from_genes.txt
	bedtools window -a ${filedir}WB_gene_positions.gff -b ${dir}/${sample}/final_results/${sample}_telocate_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\treference\t"$0}' > references_in_genes.txt

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

	#FIVE_UTR FILE
	bedtools window -a ${filedir}WB_fiveUTR_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_insertion_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tinsertion\t"$0}' > insertions_into_fiveUTR.txt
	bedtools window -a ${filedir}WB_fiveUTR_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_absence_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tabsence\t"$0}' > absences_from_fiveUTR.txt
	bedtools window -a ${filedir}WB_fiveUTR_positions.gff -b ${dir}/${sample}/final_results/${sample}_telocate_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\treference\t"$0}' > references_in_fiveUTR.txt

	#FTHREE_UTR FILE
	bedtools window -a ${filedir}WB_threeUTR_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_insertion_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tinsertion\t"$0}' > insertions_into_threeUTR.txt
	bedtools window -a ${filedir}WB_threeUTR_positions.gff -b ${dir}/${sample}/final_results/${sample}_temp_absence_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\tabsence\t"$0}' > absences_from_threeUTR.txt
	bedtools window -a ${filedir}WB_threeUTR_positions.gff -b ${dir}/${sample}/final_results/${sample}_telocate_nonredundant.bed -w 1 >intermediate.txt
	cat intermediate.txt | awk -v sample="$sample" '{print sample"\treference\t"$0}' > references_in_threeUTR.txt


	#put together
	echo "6"
	cat insertions_into_genes.txt absences_from_genes.txt references_in_genes.txt | cut -f1-3,5-7,9,11,13,15-17 >> all_gene_insertions_and_excisions.txt
	cat insertions_into_introns.txt absences_from_introns.txt references_in_introns.txt | cut -f1-3,5-7,9,11,13,15-17 >> all_introns_insertions_and_excisions.txt
	cat insertions_into_exons.txt absences_from_exons.txt references_in_exons.txt | cut -f1-3,5-7,9,11,13,15-17 >> all_exons_insertions_and_excisions.txt
	cat insertions_into_fiveUTR.txt absences_from_fiveUTR.txt references_in_fiveUTR.txt | cut -f1-3,5-7,9,11,13,15-17 >> all_fiveUTR_insertions_and_excisions.txt
	cat insertions_into_threeUTR.txt absences_from_threeUTR.txt references_in_threeUTR.txt | cut -f1-3,5-7,9,11,13,15-17 >> all_threeUTR_insertions_and_excisions.txt


done < "${file}"


python $process_script all_gene_insertions_and_excisions.txt
python $process_script all_introns_insertions_and_excisions.txt
python $process_script all_exons_insertions_and_excisions.txt
python $process_script all_fiveUTR_insertions_and_excisions.txt
python $process_script all_threeUTR_insertions_and_excisions.txt

python $check_script


cat final_all_gene_insertions_and_excisions.txt\
  final_all_introns_insertions_and_excisions.txt\
   final_all_exons_insertions_and_excisions.txt\
   final_all_fiveUTR_insertions_and_excisions.txt\
   final_all_threeUTR_insertions_and_excisions.txt > final_gene_table.txt


cat final_gene_table.txt | sort -k2,2 -k4,4 -k3,3 -k5,5n| awk '!x[$2,$3,$4,$5,$6]++' > final_gene_table_dedup.txt

#separate genes from all other events

 cat final_gene_table_dedup.txt | awk '$4=="gene" {print $0}' > dedup_genes.txt
 cat final_gene_table_dedup.txt | awk '$4!="gene" {print $0}' > dedup_nongene.txt

 python $other_gene