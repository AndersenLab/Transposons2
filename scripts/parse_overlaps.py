#!/usr/bin/env python
import re
import os
import sys
from collections import defaultdict

phenotypes="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/rnai_phenotypes.WS245.wb"
results_dir="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results"
phenos={}
phenos=defaultdict(list) # need to be able to append multiple phenos for some genes
pattern="arrest|lethal|sterile|slow growth|sick|severe"

#match phentype to transcript name
with open(phenotypes, 'r')as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('[\t]',line)
		gene_name, transcript_name, pheno=items[0:3]
		for match in re.findall(pattern, pheno):
			phenos[transcript_name].append(pheno)

#match biotype to transcript name
gene_transcripts={}
alias_codes={}
original_gene_giff="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/gene.gff"
with open(original_gene_giff, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('[\t]', line)
		WB=items[1]
		if WB=="WormBase":
			chromosome,WB,Gene_type,WB_start,WB_end,NA1,WB_orient,NA2,transcript=items[0:len(items)]
			match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", transcript) #jsut pull gene name, remove splice info
			transcript_name=match.group(1)	
			match2=re.search("biotype=(protein_coding|pseudogene|transposon_pseudogene);", transcript) #only take these 3 biotypes
			match3=re.search("Alias=(.*)",transcript)
			if match2:
				biotype=match2.group(1)
				gene_transcripts[transcript_name]=biotype
			if match3:
				aliases=re.split(",",match3.group(1))
				for alias in aliases:
					if alias != transcript_name:
						alias_codes[transcript_name]=alias


#go through overlap file
found={}
overall={}
OUT=open("essentiality_redundant.txt",'w')
with open("overlaps.txt", "r") as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('[\t]',line)
		items[0:len(items)]
		chromosome,WB,Gene_type,WB_start,WB_end,WB_orient,transcript,TE_start,TE,method,TE_class=items[0:len(items)]	
		match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", transcript) #jsut pull gene name, remove splice info
		transcript_name=match.group(1)

		ID=TE+"_"+method
		full_ID=TE+"_"+method+"_"+transcript_name #use transcript name in ID to allow for nested genes
		found[full_ID]=line
		overall[ID]=line

		if transcript_name in phenos.keys(): # add in phenotype information
			phenotype=phenos[transcript_name]
			phenotype=','.join(phenotype)
		else:
			phenotype="NA"
		match2=re.search(".*_\d+_(.*)_(\w+-)?reference",TE)
		family=match2.group(1)

		if transcript_name in gene_transcripts.keys(): # this will avoid transcript that are ncRNAs, piRNAs, excetra....manually checked transcripts failing this step to confirm
			biotype=gene_transcripts[transcript_name]
		else:
			biotype="other"
			Gene_type="gene"

		if transcript_name in alias_codes.keys():
			alias=alias_codes[transcript_name]
		else:
			alias="NA"


		OUT.write("{chromosome}\t{TE_start}\t{method}\t{family}\t{Gene_type}\t{transcript_name}\t{alias}\t{biotype}\t{phenotype}\n".format(**locals()))
		
with open("gene_overlap.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('[\t]', line)
		chromosome,WB,Gene_type,WB_start,WB_end,WB_orient,transcript,TE_start,TE,method,TE_class=items[0:len(items)]
		match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", transcript) #jsut pull gene name, remove splice info
		transcript_name=match.group(1)
		match2=re.search(".*_\d+_(.*)_(\w+-)?reference",TE)
		family=match2.group(1)
		match3=re.search("Alias=(.*)",transcript)
		if match3:
			aliases=re.split(",",match3.group(1))
			for alias in aliases:
				if alias != transcript_name:
					alias_codes[transcript_name]=alias
		ID=TE+"_"+method
		full_ID=TE+"_"+method+"_"+transcript_name
		biotype=gene_transcripts[transcript_name]
		if full_ID not in found.keys(): # if in a gene but not in exon, intron, UTR, promoter
			OUT.write("{chromosome}\t{TE_start}\t{method}\t{family}\t{Gene_type}\t{transcript_name}\t{alias}\t{biotype}\t{phenotype}\n".format(**locals()))
			overall[ID]=line

#print out TEs in intergenic regions 
with open("{results_dir}/CtCp_clipped.txt".format(**locals())) as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('[\t]', line)
		chromosome,TE_start,TE_end,TE,blank,method,strain,TE_class=items[0:len(items)]
		match2=re.search(".*_\d+_(.*)_(\w+-)?reference",TE)
		family=match2.group(1)
		ID=TE+"_"+method
		
		if ID not in overall.keys():
			transcript_name="NA"
			phenotype="NA"
			OUT.write("{chromosome}\t{TE_start}\t{method}\t{family}\tintergenic\t{transcript_name}\tNA\tNA\t{phenotype}\n".format(**locals()))
OUT.close()
#os.system("cat test.txt|sort|uniq >tmp && mv tmp test.txt")

# base call on priority order exon > five prime UTR > three prime UTR > promoter > intron > gene > intergenic
exon={}
five_prime_UTR={}
three_prime_UTR={}
promoter={}
intron={}
gene={}
intergenic={}
IDs={}
with open("essentiality_redundant.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('[\t]', line)
		chromosome,TE_start,method,family,gene_type,transcript_name,alias,biotype,phenotype=items[0:len(items)]
		ID="{chromosome}_{TE_start}_{method}_{family}_{transcript_name}".format(**locals())
		#redundant dictionaries
		if gene_type=="exon":
			exon[ID]=line
		if gene_type=="five_prime_UTR":
			five_prime_UTR[ID]=line
		if gene_type=="three_prime_UTR":
			three_prime_UTR[ID]=line
		if gene_type=="promoter":
			promoter[ID]=line
		if gene_type=="intron":
			intron[ID]=line
		if gene_type=="gene":
			gene[ID]=line
		if gene_type=="intergenic":
			intergenic[ID]=line
		if ID not in IDs.keys(): # add all IDs to a dictionary
			IDs[ID]=line

OUT=open("essentiality_nonredundant.txt", 'w')
#nonredundant classifications
for ID in IDs.keys():
	if ID in exon.keys():
		OUT.write(exon[ID]+'\n')
	elif ID in five_prime_UTR.keys():
		OUT.write(five_prime_UTR[ID]+'\n')
	elif ID in three_prime_UTR.keys():
		OUT.write(three_prime_UTR[ID]+'\n')
	elif ID in promoter.keys():
		OUT.write(promoter[ID]+'\n')
	elif ID in intron.keys():
		OUT.write(intron[ID]+'\n')
	elif ID in gene.keys():
		OUT.write(gene[ID]+'\n')
	elif ID in intergenic.keys():
		OUT.write(intergenic[ID]+'\n')
	else:
		print "ERROR, ID not classified, exiting..."
		sys.exit()
OUT.close()
#os.system("cat test2.txt| sort -k1,1 -k2,2n -k3,3 -k4,4 > tmp && mv tmp test2.txt")
