#!/bin/bash
# this script runs the kin matrix script and checks if the unique TEs are in genic or intergenic regions
# USE: GENE.sh in final_results dir

file_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files
scripts_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts
results_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/results

# pull only Wormbase protein_coding and psedogenes for gene gff file:
cat ${file_dir}/gene.gff | awk '$2=="WormBase" && $9 ~/biotype=(protein_coding|pseudogene|transposon_pseudogene)/  {print $0}'> ${file_dir}/WB_gene_positions.gff #avoid RNAs
cat  ${file_dir}/intron.gff|awk '$2=="WormBase" {print $0}' > ${file_dir}/WB_intron_positions.gff
cat  ${file_dir}/exon.gff |awk '$2=="WormBase"  {print $0}' > ${file_dir}/WB_exon_positions.gff
cat  ${file_dir}/five_prime_UTR.gff |awk '$2=="WormBase" {print $0}' > ${file_dir}/WB_fiveUTR_positions.gff
cat  ${file_dir}/three_prime_UTR.gff |awk '$2=="WormBase" {print $0}' > ${file_dir}/WB_threeUTR_positions.gff

#define promoter regions:
python ${scripts_dir}/generate_promoter_gff.py 

# run bedtools with window set to 0 for absences and references
cd ${results_dir}/final_results
mkdir gene_interrupt
cd gene_interrupt

bedtools window -a ${file_dir}/WB_gene_positions.gff -b ${results_dir}/final_results/CtCp_clipped.gff -w 0 |awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9"\t"$13"\t"$12"\t"$15"\t"$18}' > gene_overlap.txt ##the output files
bedtools window -a ${file_dir}/WB_gene_positions.gff -b ${results_dir}/final_results/CtCp_clipped.gff -w 0 |awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9"\t"$13"\t"$12"\t"$15"\t"$18}' > gene_overlap.txt ##the output files
bedtools window -a ${file_dir}/WB_exon_positions.gff -b ${results_dir}/final_results/CtCp_clipped.gff -w 0 |awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9"\t"$13"\t"$12"\t"$15"\t"$18}'  > exon_overlap.txt
bedtools window -a ${file_dir}/WB_threeUTR_positions.gff -b ${results_dir}/final_results/CtCp_clipped.gff -w 0 |awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9"\t"$13"\t"$12"\t"$15"\t"$18}'  > threeUTR_overlap.txt
bedtools window -a ${file_dir}/WB_fiveUTR_positions.gff -b ${results_dir}/final_results/CtCp_clipped.gff -w 0 |awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9"\t"$13"\t"$12"\t"$15"\t"$18}' > fiveUTR_overlap.txt
bedtools window -a ${file_dir}/WB_promoter_positions.gff -b ${results_dir}/final_results/CtCp_clipped.gff -w 0 |awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9"\t"$13"\t"$12"\t"$15"\t"$18}' > promoter_overlap.txt


#merge
cat intron_overlap.txt exon_overlap.txt threeUTR_overlap.txt fiveUTR_overlap.txt promoter_overlap.txt>  overlaps.txt

python ${scripts_dir}/parse_overlaps.py
python ${scripts_dir}/add_GO_terms.py ${file_dir}/WB_GO.txt essentiality_nonredundant.txt

cat essentiality_nonredundant_GO.txt|sort -k1,1 -k2,2n -k3,3 -k4,4 -k6,6 > tmp && mv tmp essentiality_nonredundant_GO.txt
echo  -e  "Chromosome\tTE_start\tMethod\tTE\tRegion\tTranscript_Name\tBiotype\tPhenotype\tGO_Annotation" | cat  - essentiality_nonredundant_GO.txt > tmp && mv tmp essentiality_nonredundant_GO.txt


