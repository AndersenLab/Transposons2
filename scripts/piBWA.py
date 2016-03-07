#!/usr/bin/env python
# this script: pulls the seqs of query piRNA regions, bwa aligns the seqs to the TE consensus seqeunces, and determines if the families match,
# checks whether the known Tc3 piRNAs align to the Tc3 TE sequences
# uSE: python piBWA.py (in piRNA directory)

import re
from collections import defaultdict
from subprocess import Popen, PIPE
from Bio import SeqIO
import os
import fnmatch
import sys

pairs=defaultdict(list)
found=defaultdict(list)

infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA/vc_PI.txt"
pirs="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_piRNA_positions.gff"
reference="/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa"
TE_consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta"
family_renames="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt"
known_Tc3_pi="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/known_Tc3_piRNAs.fasta"

pi_transcripts=list()
CIs=list()
# pull transcript names of the piRNA of interest
with open(infile, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		trait=items[15]
		transcript=items[22]
		pairs[trait].append(transcript)
		pi_transcripts.append(transcript)

		print trait
		CI_L=items[17]
		CI_R=items[18]
		print CI_L
		print CI_R
		#CIs[trait]
sys.exit()
# dictionary of information for piRNA transcript chrom, start, and end
with open(pirs, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		chromosome=items[0]
		start=str(int(items[3])-1) #convert to 0 base
		end=str(items[4]) #convert to 0 base
		info=items[8]
		match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", info) #jsut pull gene name, remove splice info
		transcript_name=match.group(1)	
		if transcript_name in pi_transcripts:
			found[transcript_name].extend([chromosome,start,end])
# check
for i in pi_transcripts:
	if i not in found.keys():
		print "Transcript not found, exiting"
		sys.exit()

# prep file for bedtools
OUT=open("search_regions.txt", 'w')
for i in pairs.keys():
	trait=i
	transcripts=pairs[i]
	for transcript in transcripts:
		#values=found[transcript]
		values='\t'.join(found[transcript])
		OUT.write(values + '\t' + transcript + ':' + trait + '\n')
OUT.close()


# run bedtools getfasta
cmd="bedtools getfasta -name -fi {reference} -bed search_regions.txt -fo found_piRNAs.txt".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()



#dictionary of rename list
seen=list()
renames=defaultdict(list)
with open(family_renames, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		element,family = items[0:2]
		renames[family].append(element)
		#seen.append(element)
		#seen.append(family)


#fasta_sequences = SeqIO.parse(open(TE_consensus),'fasta')
#for fasta in fasta_sequences:
#	name, sequence = fasta.id, str(fasta.seq) 
#	elements=renames[i]
#	if name not in seen:
#		seen.append(name)


# put piRNA seqs in separate files based on trait
for i in set(pairs.keys()):

	match = re.search(".*_TRANS_(.*)", i)
	transposon=match.group(1)

	OUT=open("{i}_piRNAs.fasta".format(**locals()), 'w')
	fasta_sequences = SeqIO.parse(open("found_piRNAs.txt"),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq) 
		match = re.search(":(.*_TRANS_(.*))", name) #just pull gene name, remove splice info
		trait=match.group(1)
		TE=match.group(2)
		if trait == i:
			OUT.write(">" + name + '\n' + sequence + '\n')
	OUT.close()

	# pull TE seqs of that trait
	OUT=open("{i}_TE_seqs.fasta".format(**locals()), 'w')
	fasta_sequences = SeqIO.parse(open(TE_consensus),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq) 
		if transposon in renames.keys():
			elements=renames[transposon]
		else:
			 elements="NA"
		if name == transposon or name in elements:
			OUT.write(">" + name + '\n' + sequence + '\n')
	OUT.close()

	# create bwa index for TE seq sof that family
	cmd="bwa index {i}_TE_seqs.fasta".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	# run bwa aln
	cmd= "bwa aln -o 0 -n 8 -t 2 {i}_TE_seqs.fasta {i}_piRNAs.fasta > {i}.sai".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd= "bwa samse {i}_TE_seqs.fasta {i}.sai {i}_piRNAs.fasta > {i}.sam".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd= "samtools view -bS {i}.sam > {i}.bam".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd= "samtools flagstat {i}.bam >{i}_stats.txt".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


# Check if Tc3 QTL file exists
if os.path.isfile('ZERO_new_TRANS_Tc3_TE_seqs.fasta'):
	print "Tc3 QTL file exists"
	# Test aligning the known Tc3 piRNAs to the Tc3 TE sequences
	cmd="""bwa aln -o 0 -n 2 -t 1 ZERO_new_TRANS_Tc3_TE_seqs.fasta {known_Tc3_pi} > known_Tc3.sai;
	bwa samse ZERO_new_TRANS_Tc3_TE_seqs.fasta known_Tc3.sai {known_Tc3_pi} > known_Tc3.sam;
	samtools view -bS known_Tc3.sam > known_Tc3.bam;
	samtools flagstat known_Tc3.bam > known_Tc3_stats.txt""".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	# Test aligning the known Tc3 piRNAs to the whole genome
	cmd="""bwa index {reference};
	bwa aln -o 0 -n 2 -t 1 {reference} {known_Tc3_pi} > known_Tc3_genome.sai;#
	bwa samse {reference} known_Tc3_genome.sai {known_Tc3_pi} > known_Tc3_genome.sam;#
	samtools view -bS known_Tc3_genome.sam > known_Tc3_genome.bam;
	samtools flagstat known_Tc3_genome.bam > known_Tc3_stats_genome.txt""".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()
else:
	print "Tc3 QTL file does not exist"

