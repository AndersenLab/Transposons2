#!/usr/bin/env python
# this script takes 3-4 letter gene names as input and pulls out the line of info from essentiality_nonredundant_strain_info.txt if
# the gene names matches those that were supplied
# also outputs fasta seunces on UTRs, exons, and introns 
# USE: pull_strains_info_ins.py <gene names>
# example: python ../../../scripts/pull_strains_ins_info.py prg-1


import re
import sys
from collections import defaultdict
from subprocess import Popen, PIPE
from Bio import SeqIO

queries=sys.argv[1:]
in_file="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/essentiality_nonredundant_strain_info.txt"
exon_overlap="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/exon_overlap.txt"
OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/genes_of_interest.txt", 'w')

five="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_fiveUTR_positions.gff"
three="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_threeUTR_positions.gff"
intron="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_intron_positions.gff"
exon="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_exon_positions.gff"


exons=defaultdict(list)

with open(in_file, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		chromosome=items[0]
		start=items[1]
		te=items[3]
		region=items[4]
		transcript_name=items[5]
		gene_name=items[6]
		for query in queries:
			if gene_name == query:
				OUT.write(line)
				if region == "exon":
					if gene_name in exons.keys():
						print "check for duplicates"
						sys.exit()
					else:
						exons[transcript_name].extend([chromosome,start,te,gene_name])
						#exons[gene_name].extend([chromosome,start,te,transcript_name])
OUT.close()
for i in exons.keys():
	print i

splices={}
splices_location={}
splices_gene={}
splices_transposon={}
with open(exon_overlap, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		chromosome_e,wb,region,start_e,end_e,orient,parent,te_location,te_info=items[0:9]
		match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", parent) 
		match2 = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)(.*)", parent) 
		parent_transcript=match.group(1)
		splice_form=match2.group(1)

		match = re.search("([A-Z]+)_(\d+)_(([A-Za-z\d+_-]+))_((\w+-)?reference)", te_info) #([A-Za-z\d+_])_((\w+-)?reference)\w+_\d+_\d+
		test=match.group(1)
		transposon_location=match.group(2)
		transposon=match.group(3)

		if parent_transcript in exons.keys():

			info=exons[parent_transcript]
			C1,S1,T1,G1=info[0:len(info)]
			if chromosome_e == C1: # check that chromosome matches
				if transposon_location==S1: # check that location of insertion matches
					if transposon == T1: # check that transposon family matches


						splices[splice_form]=parent_transcript
						splices_location[splice_form]=str(te_location)
						splices_gene[splice_form]=G1
						splices_transposon[splice_form]=transposon

# get lengths of genes based on parent transcripts
gene_lengths={}
gene_starts={}
gene_gff="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_gene_positions.gff"
with open(gene_gff, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		chromsome, wb,region,start,end,blank,orient,blank2,parent=items[0:9]
		match2 = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)(.*);biotype", parent) 
		gene=match2.group(1)
		if gene in splices.values():
			length=int(end)-int(start)
			gene_lengths[gene]=length
			gene_starts[gene]=start


five_splice=defaultdict(list)
three_splice=defaultdict(list)
exon_splice=defaultdict(list)
intron_splice=defaultdict(list)

# five prime
with open(five, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		chromsome, wb,region,start,end,blank,orient,blank2,parent=items[0:9]
		match2 = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)(.*)", parent) 
		splice_form=match2.group(1)
		if splice_form in splices.keys():
			five_splice[splice_form].extend([chromsome,start,end])

# three prime
with open(three, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		chromsome, wb,region,start,end,blank,orient,blank2,parent=items[0:9]
		match2 = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)(.*)", parent) 
		splice_form=match2.group(1)
		if splice_form in splices.keys():
			three_splice[splice_form].extend([chromsome,start,end])			

# exon 
with open(exon, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		chromsome, wb,region,start,end,blank,orient,blank2,parent=items[0:9]
		match2 = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)(.*)", parent) 
		splice_form=match2.group(1)
		if splice_form in splices.keys():
			exon_splice[splice_form].extend([chromsome,start,end])	

# intron 
with open(intron, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		chromsome, wb,region,start,end,blank,orient,blank2,parent=items[0:9]
		match2 = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([A-Za-z\.\d+_-]+);?", parent) 
		splice_form=match2.group(1)
		if splice_form in splices.keys():
			intron_splice[splice_form].extend([chromsome,start,end])	


#write out region info to a new file	
for i in splices.keys():
	outfile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/models/gene_model_regions{i}.txt".format(**locals())
	fasta="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/models/gene_model_regions{i}.fasta".format(**locals())
	OUT=open(outfile, 'w')
	five_regions=five_splice[i]
	for x in range(0,len(five_regions),3):
		xx=x/3
		OUT.write(five_regions[x] + '\t')
		OUT.write(five_regions[x+1] + '\t')
		OUT.write(five_regions[x+2] + "\tFIVE_{i}_{xx}\n".format(**locals()))
		#OUT.write(five_regions[x+2] + '\t' + i +'_' + str(x) +' \n')
	exon_regions=exon_splice[i]
	for x in range(0,len(exon_regions),3):
		xx=x/3
		OUT.write(exon_regions[x] + '\t')
		OUT.write(exon_regions[x+1] + '\t')
		OUT.write(exon_regions[x+2] + "\tEXON_{i}_{xx}\n".format(**locals()))
	intron_regions=intron_splice[i]
	for x in range(0,len(intron_regions),3):
		xx=x/3
		OUT.write(intron_regions[x] + '\t')
		OUT.write(intron_regions[x+1] + '\t')
		OUT.write(intron_regions[x+2] + "\tINTRON_{i}_{xx}\n".format(**locals()))
	three_regions=three_splice[i]
	for x in range(0,len(three_regions),3):
		xx=x/3
		OUT.write(three_regions[x] + '\t')
		OUT.write(three_regions[x+1] + '\t')
		OUT.write(three_regions[x+2] + "\tTHREE_{i}_{xx}\n".format(**locals()))	
	OUT.close()

	#sort the bed file
	cmd="bedtools sort -i {outfile} > tmp.txt && mv tmp.txt {outfile}".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


	#run bedtools getfasta
	cmd="bedtools getfasta -name -fi /lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa -bed {outfile} -fo {fasta}".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/models/gene_model_regions{i}.model".format(**locals()), 'w')

	location=splices_location[i]
	parent=splices[i]
	parent_start=gene_starts[parent]
	marker=int(location)-int(parent_start)



	OUT.write(i +"\t" + str(marker) + "\t" + splices_gene[i] + "\t" + splices_transposon[i] + "\n")
	total_length=0
	#print exons in caps, introns in lower case, and UTRs on their own lines
	fasta_sequences = SeqIO.parse(open(fasta),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq) 
		match=re.search("(\w+)_.*", name)
		region=match.group(1)
		if region=="EXON":
			OUT.write(sequence.upper())
			total_length +=len(sequence)
		elif region=="INTRON":
			OUT.write(sequence)
			total_length +=len(sequence)
		else:
			OUT.write('\n>')
			OUT.write(name)
			OUT.write('\n')
			OUT.write(sequence)
			OUT.write('\n')
	OUT.close()

	

