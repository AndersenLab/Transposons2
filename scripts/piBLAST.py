#!/usr/bin/env python
# this script pull the seqs of query piRNA regions, blasts the seqs to the TE consensus seqeunces, and determines if the families match
# uSE: python piBLAST.py (in piRNA directory)

import re
from collections import defaultdict
from subprocess import Popen, PIPE

pairs={}
found=defaultdict(list)

infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA/vc_PI.txt"
pirs="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_piRNA_positions.gff"
reference="/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa"
TE_consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta"
family_renames="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt"


# pull transcript names of the piRNA of interest
with open(infile, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		
		trait=items[15]
		transcript=items[22]
		pairs[transcript]=trait

# dictionary of information for piRNA transcript chrom, start, and end
with open(pirs, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		chromosome=items[0]
		start=str(int(items[3])-1) #convert to 0 base
		end=str(items[4])
		info=items[8]
		match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", info) #jsut pull gene name, remove splice info
		transcript_name=match.group(1)	
		if transcript_name in pairs.keys():
			found[transcript_name].extend([chromosome,start,end])
# check
for i in pairs.keys():
	if i not in found.keys():
		print "Transcript not found, exiting"
		sys.exit()

# prep file for bedtools
OUT=open("search_regions.txt", 'w')
for i in found.keys():
	trait=pairs[i]
	values='\t'.join(found[i])
	OUT.write(values + '\t' + i + ':' + trait + '\n')
OUT.close()


# run bedtools getfasta
cmd="bedtools getfasta -name -fi {reference} -bed search_regions.txt -fo found_piRNAs.txt".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

# make blast database of TE sequences
cmd="/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/makeblastdb -in  {TE_consensus} -dbtype nucl -out TE_database".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

# blast piRNA seqeunces to TE sequences
cmd="/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/blastn -db TE_database -query found_piRNAs.txt -evalue 1 -word_size 5 -outfmt 6 -max_target_seqs 100 -out piRNA_blast.txt -num_threads 10"
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

# correct redundant names/merge TEs
cmd="python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/correct_names_2_set2.py {family_renames} piRNA_blast.txt TMP".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()
cmd="mv TMP piRNA_blast.txt".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


# check if the trait family matched the family hit by the piRNA blast
OUT=open("blast_results.txt", 'w')
with open("piRNA_blast.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		result_fam=items[1]
		match = re.search("TRANS_(.*)", items[0]) #just pull gene name, remove splice info
		TE=match.group(1)	

		if TE == result_fam:
			decision="MATCHED"
		else:
			decision="--"
		OUT.write(line + '\t' + decision + '\n')
OUT.close()








