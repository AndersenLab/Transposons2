#!/usr/bin/env python
# this scrip check if any piRNA sequences align to the TE sequences using various numbers of allowed mismatches
# USE: pi_v_TE.py

import re
import sys
import pysam
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


OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/piRNA_regions.bed", 'w')
with open(pi_IN, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		chromsome,WB,pi,start,end=items[0:5]
		orient=items[6]
		piRNA=items[8]
		start=int(start)-1
		OUT.write("{chromsome}\t{start}\t{end}\t{piRNA}\t.\t{orient}\n".format(**locals()))
OUT.close()
# get fasta seqs of piRNAs
cmd="bedtools getfasta -s -name -fi {reference} -bed /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/piRNA_regions.bed -fo {pi_fasta}".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

# create bwa index for TE seqs of that family
cmd="bwa index {TE_consensus}".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

OUT_SUMMARY=open("summary_mismatches_BWA.txt", 'w')
OUT_SUMMARY.write("Number of Mismatches\tNumber Unique piRNAs Aligned\tNumber Unique Transposons\n")

OUT_SUMMARY_STRICT=open("summary_mismatches_BWA_strict.txt", 'w')
OUT_SUMMARY_STRICT.write("Number of Mismatches\tNumber Unique piRNAs Aligned\tNumber Unique Transposons\n")

def align(mismatches,strict=False):
	TE_consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta"
	pi_fasta="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/piRNAs.fasta"
	family_renames="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt"
	
	# run bwa aln
	cmd = "bwa aln -o 0 -n {mismatches} -t 2 {TE_consensus} {pi_fasta} > {mismatches}_pi_v_TE.sai".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd = "bwa samse {TE_consensus} {mismatches}_pi_v_TE.sai {pi_fasta} > {mismatches}_pi_v_TE.sam".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


	# rename headers and read info
	OUT=open("{mismatches}_renamed.sam".format(**locals()), 'w')
	with open("{mismatches}_pi_v_TE.sam".format(**locals()) ,'r') as IN:
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


	cmd = "samtools view -bS -F4 {mismatches}_renamed.sam > {mismatches}_renamed.bam".format(**locals()) # filter out unmapped reads
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd = "samtools sort -o -@ 8 {mismatches}_renamed.bam out > {mismatches}_renamed.sorted.bam".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd = "samtools index {mismatches}_renamed.sorted.bam".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate() 

	cmd = "samtools flagstat {mismatches}_renamed.sorted.bam > {mismatches}_stats.txt".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd = "samtools view {mismatches}_renamed.sorted.bam |cut -f1|sort|uniq |wc -l".format(**locals())
	unique_pis, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	cmd = "samtools view {mismatches}_renamed.sorted.bam |cut -f3|sort|uniq |wc -l".format(**locals())
	unique_TEs, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

	unique_pis=re.sub('\n','',unique_pis)
	unique_TEs=re.sub('\n','',unique_TEs)

	OUT_SUMMARY.write("{mismatches}\t{unique_pis}\t{unique_TEs}\n".format(**locals()))


	Bfile = pysam.AlignmentFile("{mismatches}_renamed.sorted.bam".format(**locals()), "rb")

	seen_pis={}
	seen_TEs={}
	seen_pis_strict={}
	seen_TEs_strict={}

	TEST=open("test.txt",'w')
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
		OUT_SUMMARY_STRICT.write("{mismatches}\t{no_pis_strict}\t{no_TEs_strict}\n".format(**locals()))

	no_pis = len(seen_pis.keys())
	no_TEs = len(seen_TEs.keys())


	# make sure methods of counting uniqueness are the same
	if int(no_pis) != int(unique_pis):
		sys.exit("ERROR: Inconsistency in unique counts,exiting...")
	if int(no_TEs) != int(unique_TEs):
		sys.exit("ERROR: Inconsistency in unique counts,exiting...")


align(0,strict=True)
align(1,strict=True)
align(2,strict=True)
align(3,strict=True)

OUT_SUMMARY.close()
OUT_SUMMARY_STRICT.close()


