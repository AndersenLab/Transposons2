#!/usr/bin/env python
# this script: pulls the seqs of query piRNA regions, bwa aligns the seqs to the TE consensus seqeunces, and determines if the families match,
# checks whether the known Tc3 piRNAs align to the Tc3 TE sequences
# check if any trait with a QTL in a piRNA regions has piRNA within its CI that align to the corresponding TE seqs
# USE: python piBWA.py (in piRNA directory)
# NOTE: must first transfer "vc_PI.txt"

import re
from collections import defaultdict
from subprocess import Popen, PIPE
from Bio import SeqIO
import os
import fnmatch
import sys
import pickle
import pysam

pairs=defaultdict(list)
found=defaultdict(list)
CIs=defaultdict(list)
pairs_CI={}
found_CIs=defaultdict(list)
pi_transcripts=list()

infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA/vc_PI.txt"
piqtl="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/QLT_piRNA_overlap.txt" #"/lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA/TT.txt"
pirs="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_piRNA_positions.gff"
reference="/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa"
TE_consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta"
family_renames="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt"
known_Tc3_pi="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/known_Tc3_piRNAs.fasta"


shorts={}
# pull transcript names of the piRNA of interest
with open(infile, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		chromosome=items[0]
		trait=items[14]
		transcript=items[21]
		pairs[trait].append(transcript)
		pi_transcripts.append(transcript)

# put shortened WB family names into a dictionary
with open(family_renames, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		family=items[1]
		family_short=re.sub("_CE$","",family)
		family_short=re.sub("WBTransposon","WBT",family_short)
		shorts[family_short]=family


# put shortened consensus family names into a dictionary
fasta_sequences = SeqIO.parse(open(TE_consensus),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq) 
	family_short=re.sub("_CE$","",name)
	family_short=re.sub("WBTransposon","WBT",family_short)
	shorts[family_short]=name



# CI version, create bed file
OUT=open("CIs.bed", 'w')
with open(piqtl, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		
		trait=items[0]
		chromosome=items[2]
		CI_L=items[6]
		CI_R=items[7]


		
		trait=re.sub('(ins)', 'ONE_new',trait)
		trait=re.sub('(ref)', 'reference',trait)
		trait=re.sub('(abs)', 'absent',trait)

		match=re.search('(.*)\((.*)\)',trait)
		fam=match.group(1)
		method=match.group(2)

		if fam in shorts.keys():
			fam=shorts[fam]

		trait=method + "_TRANS_" + fam
		pairs_CI[trait]=0
		OUT.write("{chromosome}\t{CI_L}\t{CI_R}\t{trait}\t.\t.\n".format(**locals()))
OUT.close()
print "STEP 2"
# find interesection of CIs and piRNAs
if not os.path.isfile("CIs_pirs_intersect.txt"):
	# find interesection of CIs and piRNAs
	cmd = "bedtools intersect -wo -a CIs.bed -b {pirs} > CIs_pirs_intersect.txt".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


	pi_position_info=defaultdict(list)
	pi_trait_info=defaultdict(list)
	#iterate over intersection file and write search region file for the CIs to prep for bedtools
	#iterate over intersection file and write search region file for the CIs to prep for bedtools
	
	with open("CIs_pirs_intersect.txt", 'r') as IN:
		for line in IN:
			line=line.rstrip('\n')
			items=re.split('\t', line)
			trait_hit=items[3]
			chromosome, start, end, orient = items[6],items[9],items[10],items[12]
			info=items[14]
			match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", info) #just pull gene name, remove splice info
			transcript_name=match.group(1)	
			pi_position_info[transcript_name]=[chromosome,start,end,orient]
			pi_trait_info[transcript_name].append(trait_hit)


	with open("pi_position_info.txt", "wb") as fp: # Pickle
		pickle.dump(pi_position_info, fp) 

	with open("pi_trait_info.txt", "wb") as fp: # Pickle
		pickle.dump(pi_trait_info, fp) 



else:
	print "CIs_pirs_intersect.txt exists, continuing..."
	with open("pi_position_info.txt", "rb") as fp: 
		pi_position_info = pickle.load(fp)  # Unpickle
	with open("pi_trait_info.txt", "rb") as fp: 
		pi_trait_info = pickle.load(fp)  # Unpickle



print "STEP 3"
OUT=open("search_regions_fullPI.txt",'w')
for i,z in pi_position_info.items():
	chromosome, start, end, orient = z
	start=int(start)-1 #convert to 0 base
	traits = set(pi_trait_info[i])
	OUT.write("{chromosome}\t{start}\t{end}\t{i}:".format(**locals()))
	for trait in traits:
		OUT.write(trait + ";")
	OUT.write("\t.\t{orient}".format(**locals()))
	OUT.write('\n')
OUT.close()

print "STEP 4"
# dictionary of information for piRNA transcript chrom, start, and end
with open(pirs, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		chromosome=items[0]
		start=str(int(items[3])-1) #convert to 0 base
		end=str(items[4]) #convert to 0 base
		orient=items[6]
		info=items[8]
		match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", info) #jsut pull gene name, remove splice info
		transcript_name=match.group(1)	
		if transcript_name in pi_transcripts:
			found[transcript_name].extend([chromosome,start,end,orient])

print "STEP 5"
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
		values='\t'.join((found[transcript])[0:3])
		orient=(found[transcript])[3]
		OUT.write(values + '\t' + transcript + ':' + trait + '\t.\t' + orient + '\n')
OUT.close()

print "STEP 6"
# run bedtools getfasta
cmd="bedtools getfasta -name -s -fi {reference} -bed search_regions.txt -fo found_piRNAs.txt".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

# run bedtools getfasta for CIs
cmd="bedtools getfasta -name -s -fi {reference} -bed search_regions_fullPI.txt -fo found_full_piRNAs.txt".format(**locals())
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

OUT_SUMMARY=open("summary_mismatches_BWA_vi.txt", 'w')
OUT_SUMMARY.write("Number of Mismatches\tNumber Unique piRNAs Aligned\tNumber Unique Transposons\n")

OUT_SUMMARY_STRICT=open("summary_mismatches_BWA_strict_vi.txt", 'w')
OUT_SUMMARY_STRICT.write("Number of Mismatches\tNumber Unique piRNAs Aligned\tNumber Unique Transposons\n")

print "STEP 7"
def align(mismatches,pi,strict=False,full=False):
	TE_consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta"
	pi_fasta="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/piRNAs.fasta"
	family_renames="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt"
	print pi
	if full==True:
		full="_full"
	else:
		full=""

	# run bwa aln
	cmd= "bwa aln -o 0 -n {mismatches} -t 2 {pi}_TE_seqs.fasta {pi}{full}_piRNAs.fasta > {pi}{full}.sai".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd= "bwa samse {pi}{full}_TE_seqs.fasta {pi}{full}.sai {pi}{full}_piRNAs.fasta > {pi}{full}.sam".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


	# rename headers and read info
	OUT=open("{pi}{full}_renamed.sam".format(**locals()), 'w')
	with open("{pi}{full}.sam".format(**locals()) ,'r') as IN:
		for line in IN:
			line=line.rstrip()
			items=re.split('\t',line)
			TE=items[2]
			if TE in renames.keys():
				TE=renames[TE]
				items[2]=TE

			if re.search('^@SQ',line):
				sn=items[1]
				match=re.search("SN:(.*)",sn)
				element=match.group(1)
				if element in renames.keys():
					trans=renames[element]
					items[1]="SN:"+ trans
			new_line='\t'.join(items[0:])
			OUT.write(new_line + '\n')
	OUT.close()


	cmd = "samtools view -bS -F4 {pi}{full}_renamed.sam > {pi}{full}_renamed.bam".format(**locals()) # filter out unmapped reads
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd = "samtools sort -o -@ 8 {pi}{full}_renamed.bam out > {pi}{full}_renamed.sorted.bam".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd = "samtools index {pi}{full}_renamed.sorted.bam".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate() 

	cmd = "samtools flagstat {pi}{full}_renamed.sorted.bam > {pi}{full}_stats.txt".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd = "samtools view {pi}{full}_renamed.sorted.bam |cut -f1|sort|uniq |wc -l".format(**locals())
	unique_pis, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd = "samtools view {pi}{full}_renamed.sorted.bam |cut -f3|sort|uniq |wc -l".format(**locals())
	unique_TEs, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	unique_pis=re.sub('\n','',unique_pis)
	unique_TEs=re.sub('\n','',unique_TEs)

	OUT_SUMMARY.write("{pi}\t{unique_pis}\t{unique_TEs}\n".format(**locals()))






	Bfile = pysam.AlignmentFile("{pi}{full}_renamed.sorted.bam".format(**locals()), "rb")

	seen_pis={}
	seen_TEs={}
	seen_pis_strict={}
	seen_TEs_strict={}


	Binfo = Bfile.fetch()
	for x in Binfo:
		query = x.query_name
		TE = Bfile.getrname(x.reference_id)
		flag = x.flag
		MD = x.get_tag('MD')
		TEST.write(query + '\n')

		if query=="*":
			print "TTTTT"

		MD_nums=re.findall('\d+', MD)
		first_digit=int(MD_nums[0])
		last_digit=int(MD_nums[-1])

		seen_pis[query]=0
		seen_TEs[TE]=0


	# enforce no mismatches in first 8 bps of the piRNA
		if strict:
			if flag==0 and first_digit>=8:
				seen_pis_strict[query]=0
				seen_TEs_strict[TE]=0

			elif flag==16 and last_digit>=8:
				seen_pis_strict[query]=0
				seen_TEs_strict[TE]=0

			elif flag !=0 and flag !=16:
				sys.exit("ERROR: Flag  %s not accounted for, exiting..." %flag)

	if strict:
		no_pis_strict = len(seen_pis_strict.keys())
		no_TEs_strict = len(seen_TEs_strict.keys())
		OUT_SUMMARY_STRICT.write("{mismatches}\t{pi}\t{no_pis_strict}\t{no_TEs_strict}\n".format(**locals()))

	no_pis = len(seen_pis.keys())
	no_TEs = len(seen_TEs.keys())


	# make sure methods of counting uniqueness are the same
	if int(no_pis) != int(unique_pis):
		sys.exit("ERROR: Inconsistency in unique counts,exiting...")
	if int(no_TEs) != int(unique_TEs):
		sys.exit("ERROR: Inconsistency in unique counts,exiting...")





# put piRNA seqs in separate files based on trait
for i in set(pairs.keys()):
	if not re.search('total', i):
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



		# create bwa index for TE seqs of that family
		cmd="bwa index {i}_TE_seqs.fasta".format(**locals())
		result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()



		align(8,i,strict=True)













OUT_SUMMARY.close()
OUT_SUMMARY_STRICT.close()












########################
########################
# CI VERSION
OUT_SUMMARY=open("summary_mismatches_BWA_ci.txt", 'w')
OUT_SUMMARY.write("Number of Mismatches\tNumber Unique piRNAs Aligned\tNumber Unique Transposons\n")

OUT_SUMMARY_STRICT=open("summary_mismatches_BWA_strict_ci.txt", 'w')
OUT_SUMMARY_STRICT.write("Number of Mismatches\tNumber Unique piRNAs Aligned\tNumber Unique Transposons\n")
for i in set(pairs_CI.keys()):
	print i

	match = re.search(".*_TRANS_(.*)", i)
	transposon=match.group(1)

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


	OUT=open("{i}_full_piRNAs.fasta".format(**locals()), 'w')
	fasta_sequences = SeqIO.parse(open("found_full_piRNAs.txt"),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		name=name.rstrip('\n')
		#name=name.rstrip(';')
		all_traits=re.split(':', name)
		Transcript=all_traits[0]
		all_traits=all_traits[1]
		all_traits=re.split(';', all_traits)
		print name
		for x in all_traits:
			if x != '\n' and x!='':
				#x=x.rstrip('\n')
				print "TRAIT IS " + x
				match = re.search("(.*_TRANS_(.*))", x) #just pull gene name, remove splice info
				trait=match.group(1)
				if trait == i:
					OUT.write(">" + Transcript + ":" + x + '\n' + sequence + '\n')
	OUT.close()


	# create bwa index for TE seqs of that family
	cmd="bwa index {i}_TE_seqs.fasta".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	align(8,i,strict=True,full=True)

	# run bwa aln
	#cmd= "bwa aln -o 0 -n 3 -t 2 {i}_TE_seqs.fasta {i}_full_piRNAs.fasta > full_{i}.sai".format(**locals())
	#result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	#cmd= "bwa samse {i}_TE_seqs.fasta full_{i}.sai {i}_full_piRNAs.fasta > full_{i}.sam".format(**locals())
	#result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	#cmd= "samtools view -bS full_{i}.sam > full_{i}.bam".format(**locals())
	#result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	#cmd= "samtools flagstat full_{i}.bam > full_{i}_stats.txt".format(**locals())
	#result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()
OUT_SUMMARY.close()
OUT_SUMMARY_STRICT.close()



sys.exit()

########################
# Check if Tc3 QTL file exists
if os.path.isfile('ONE_new_TRANS_Tc3_TE_seqs.fasta'):
	print "Tc3 QTL file exists"
	# Test aligning the known Tc3 piRNAs to the Tc3 TE sequences
	cmd="""bwa aln -o 0 -n 3 -t 1 ONE_new_TRANS_Tc3_TE_seqs.fasta {known_Tc3_pi} > known_Tc3.sai;
	bwa samse ONE_new_TRANS_Tc3_TE_seqs.fasta known_Tc3.sai {known_Tc3_pi} > known_Tc3.sam;
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
else:
	print "Tc3 QTL file does not exist"

# generate summary file of stats
OUT=open("stats_summary.txt", 'w')
for stat_file in os.listdir('.'):
    if fnmatch.fnmatch(stat_file, '*stats.txt'):
        trait, file_extension = os.path.splitext(stat_file)
        with open(stat_file, 'r') as IN:
        	total=next(IN)
        	total_PIs=re.split('\s+', total)[0]
        	next(IN)
        	mapped=next(IN)
        	mapped_PIs=re.split('\s+', mapped)[0]
        	OUT.write(trait + '\t' + mapped_PIs + '\t' + total_PIs + '\n')
OUT.close()



########################
# Align all piRNAs to all TE seqs
#cmd="""bwa aln -o 0 -n 3 -t 1 ZERO_new_TRANS_Tc3_TE_seqs.fasta {known_Tc3_pi} > known_Tc3.sai;




