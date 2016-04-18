#!/usr/bin/env python
# this scrip check if any piRNA sequences BLAST to the TE sequences 
# USE: pi_BLAST_all.py

import re
import sys
from subprocess import Popen, PIPE 

pi_IN="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_piRNA_positions.gff"
reference="/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa"
pi_fasta="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/piRNAs.fasta"
TE_consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta"
family_renames="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt"


# put shortened WB family names into a dictionary
renames={}
with open(family_renames, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		element,family=items[0:2]
		renames[element]=family


# redundant-have already retrieved the fasta seq of piRNAs in the previosu script
OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/piRNA_regions.bed", 'w')
with open(pi_IN, 'r') as IN:
	for line in IN:
		print line
		line=line.rstrip('\n')
		items=re.split('\t',line)
		chromsome,WB,pi,start,end=items[0:5]
		orient=items[6]
		piRNA=items[8]
		start=int(start)-1
		print orient
		OUT.write("{chromsome}\t{start}\t{end}\t{piRNA}\t.\t{orient}\n".format(**locals()))
OUT.close()
# get fasta seqs of piRNAs
cmd="bedtools getfasta -s -name -fi {reference} -bed /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/piRNA_regions.bed -fo {pi_fasta}".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


# make blast database of TE sequences
cmd="/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/makeblastdb -in  {TE_consensus} -dbtype nucl -out TE_database".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

# blast piRNA seqeunces to TE sequences
cmd="/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/blastn -db TE_database -query {pi_fasta} -evalue 1 -word_size 5 -outfmt 6 -max_target_seqs 100 -out piRNA_blast.txt -num_threads 10".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


#ID=Transcript:W04C9.7;Parent=Gene:WBGene00166632;Name=W04C9.7;locus=21ur-8411	WBTransposon00000671	100.00	12	0	0	5	16	2989	3000	0.84	23.3



OUT=open("TMP.txt", 'w')
with open("piRNA_blast.txt".format(**locals()) ,'r') as IN:
	for line in IN:
		line=line.rstrip()
		items=re.split('\t',line)
		TE=items[1]
		if TE in renames.keys():
			TE=renames[TE]
			items[1]=TE
		new_line='\t'.join(items[0:])
		OUT.write(new_line + '\n')
OUT.close()

cmd="mv TMP.txt piRNA_blast.txt".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

