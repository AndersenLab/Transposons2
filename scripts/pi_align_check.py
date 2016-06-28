#!/usr/bin/env python
# this script checks if the piRNA within the CI of the QTL for a given trait aligned to the TE elements/family for that trait
# also checks if knwon Tc3 seqs from literature align to the Tc3 sequences and  to the genome
#USE: pi_align_check.py

import re
import pickle
from collections import defaultdict
from Bio import SeqIO
from subprocess import Popen, PIPE


def search_aligned(mismatch):
	OUT=open("intervalPIs_{mismatch}.txt".format(**locals()),'w') 
	qtl_overlap="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/QLT_piRNA_overlap.txt"

	in_interval=defaultdict(list)
	aligned_in_interval=defaultdict(list)

	with open("strict_alignments_3.txt", "rb") as fp: #load information on which piRNAs aligned to which TEs
		alignments = pickle.load(fp)

	with open(qtl_overlap, 'r') as IN:
		next(IN)
		for line in IN:
			line = line.rstrip('\n')
			items=re.split('\t', line)
			
			TE=items[0]
			transcript=items[8]
			
			match=re.search('(.*)\(.*\)',TE)
			family=match.group(1)

			pis_aligned=alignments[family]
			in_interval[TE].append(transcript)

			if transcript in pis_aligned:
				aligned_in_interval[TE].append(transcript)


	for i in in_interval.keys():
		no_in_interval=len(in_interval[i])
		if i in aligned_in_interval.keys():
			no_aligned_in_interval=len(aligned_in_interval[i])
		else:
			no_aligned_in_interval=0
		OUT.write("{i}\t{mismatch}\t{no_aligned_in_interval}\t{no_in_interval}\n".format(**locals()))
	OUT.close()

search_aligned(0)
search_aligned(1)
search_aligned(2)
search_aligned(3)

cmd="cat intervalPIs_* > intervalPIs_all.txt"
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()
########################
# TC3 control
known_Tc3_pi="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/known_Tc3_piRNAs.fasta"
reference="/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa"
TE_consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta"
family_renames="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt"

renames={}
# put shortened WB family names into a dictionary
with open(family_renames, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		element=items[0]
		family=items[1]
		renames[element]=family


# pull TE seqs of that trait
OUT=open("TC3_TE_seqs.fasta".format(**locals()), 'w')
fasta_sequences = SeqIO.parse(open(TE_consensus),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq) 
	if name in renames.keys():
		transposon=renames[name]
	else:
		transposon=name
	if transposon == "Tc3":
		OUT.write(">" + name + '\n' + sequence + '\n')
OUT.close()

cmd="bwa index TC3_TE_seqs.fasta".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()
# Test aligning the known Tc3 piRNAs to the Tc3 TE sequences
cmd="""bwa aln -o 0 -n 0 -t 1 TC3_TE_seqs.fasta {known_Tc3_pi} > known_Tc3.sai;
bwa samse TC3_TE_seqs.fasta known_Tc3.sai {known_Tc3_pi} > known_Tc3.sam;
samtools view -bS known_Tc3.sam > known_Tc3.bam;
samtools flagstat known_Tc3.bam > known_Tc3_stats.txt""".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

# Test aligning the known Tc3 piRNAs to the whole genome
cmd="""bwa index {reference};
bwa aln -o 0 -n 3 -t 1 {reference} {known_Tc3_pi} > known_Tc3_genome.sai;#
bwa samse {reference} known_Tc3_genome.sai {known_Tc3_pi} > known_Tc3_genome.sam;#
samtools view -bS known_Tc3_genome.sam > known_Tc3_genome.bam;
samtools flagstat known_Tc3_genome.bam > known_Tc3_genome_stats.txt""".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


########################
# check traits of interest
trait_list=[]
with open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/TOI.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		trait_list.append(line)

OUT=open("TOIs_matching_intervalPIs.txt",'w')  # file for TOIs with piRNAs found in the CI of the QTL aligning to their repective TEs  
with open("intervalPIs_all.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		trait,mismatch,found,possible=items[0:]
		found=int(found)
		if trait in trait_list:
			if found !=0:
				OUT.write(line +'\n')
OUT.close()

