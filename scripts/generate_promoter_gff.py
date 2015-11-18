#!/usr/bin/env python
# this script extracts promoter regions from the WB gff file
# promoter here is defined as the region 2000 bp upstream of the start of the gene, or up until the closest feature if that feature is within 2000 bps
# USE: generate_promoter_gff.py
import re
import sys
import os
from collections import defaultdict

file_dir="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files"
WB_seq_gff="{file_dir}/WB_gene_seq_name.gff".format(**locals())
WB_gene_gff="{file_dir}/WB_gene_positions.gff".format(**locals())

#simplify gene gff
os.system("""cat {WB_gene_gff} |awk -v OFS='\t' '{{match($9,"sequence_name=[a-zA-Z0-9._]+",a)}} {{print $1,$2,$3,$4,$5,$6,$7,$8,a[0]}}'> {WB_seq_gff}""".format(**locals()))

points="{file_dir}/WB_gene_points.gff".format(**locals()) # all starts and ends
plus_minus="{file_dir}/WB_gene_plus_minus.gff".format(**locals()) # starts of plus and ends of minus
promoters="{file_dir}/WB_promoter_positions.gff".format(**locals())


POINTS=open(points, 'w')
PLUS_MINUS=open(plus_minus, 'w')

# split up gene gff into points
with open(WB_seq_gff, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		start=items[3]
		end=items[4]
		orient=items[6]
		if orient=="+":
			PLUS_MINUS.write("\t".join(items[0:3])+"\t" + start + "\t" + start + "\t" + "\t".join(items[5:9]) + "\n") # pull the start poistions of plus oriented genes
		if orient=="-":
			PLUS_MINUS.write("\t".join(items[0:3])+"\t" + end + "\t" + end + "\t" + "\t".join(items[5:9]) + "\n") # pull the minus poistions of plus oriented genes
		POINTS.write("\t".join(items[0:3])+"\t" + start + "\t" + start + "\t" + "\t".join(items[5:9]) + "\n") # output all start positions
		POINTS.write("\t".join(items[0:3])+"\t" + end + "\t" + end + "\t" + "\t".join(items[5:9]) + "\n") # output all minus positions


POINTS.close()
PLUS_MINUS.close()

#resort points file
os.system("cat {points}  | sort -k1,1 -k4,4n > tmp && mv tmp {points}".format(**locals()))
#resort plus_minus files
os.system("cat {plus_minus}  | sort -k1,1 -k4,4n > tmp && mv tmp {plus_minus}".format(**locals()))
#run bedtools closest
os.system("bedtools closest -a {plus_minus} -b {points}  -id -D a -t all -k 3 >gene_to_gene.txt".format(**locals())) # allow 3 closeset hits to account for duplicate transcript names matching each other twice

# define promoter regions
promoter_regions=defaultdict(list)
with open("gene_to_gene.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('[\t]',line)
		chromosome=items[0]
		start1=int(items[3])
		end1=items[4]
		orient=items[6]
		seq=items[8]
		start2=items[12]
		end2=items[13]
		matched_seq=items[17]
		distance=int(items[18])
		if seq != matched_seq:
			if distance != 0: # do not want the expact overlaps (start of one gene overlapping with end of another gene) when calculating promoter regions
				if distance < -2000: # if more than 2000 bp away
					if orient == "+":
						promoter_start=start1-2000	
					if orient == "-":
						promoter_start=start1+2000

				if distance >= -2000: # if more than 2000 bp away
					if orient == "+":
						promoter_start=start1+(distance+1) # add 1 because distance is is negative....have to adjust by one so that promoter does not overlap with the gene that was closest
					if orient == "-":
						promoter_start=start1-(distance+1) # add 1 because distance is is negative....have to adjust by one so that promoter does not overlap with the gene that was closest
				if seq in promoter_regions.keys(): # if duplicate transcript name (occurs in 3 cases), take the earlier start position (only changes promoter for one gene)
					prev_start=(promoter_regions[seq])[2]
					if orient == "+":
						if start1<prev_start:
							promoter_regions[seq].extend([chromosome,orient,start1,promoter_start,distance])
					if orient == "-":
						if start1>prev_start:
							promoter_regions[seq].extend([chromosome,orient,start1,promoter_start,distance])

				else:
					promoter_regions[seq].extend([chromosome,orient,start1,promoter_start,distance])

#manually add ends of chromosomes (b/c closest hit upstream will only be to itself),note: gff 1 based:
promoter_regions["sequence_name=2L52.1"].extend(["II","+",1867,1867-1,0])
promoter_regions["sequence_name=2RSSE.3"].extend(["II","-",15278575,15278575+2000,0])
promoter_regions["sequence_name=4R79.1"].extend(["IV","-",17490497,17490497+2000,0])
promoter_regions["sequence_name=cTel3X.1"].extend(["V","+",1480,1480-1,0])
promoter_regions["sequence_name=cTel7X.1"].extend(["X","+",1316,1316-1,0])
promoter_regions["sequence_name=MTCE.3"].extend(["MtDNA","+",113,1,0])
promoter_regions["sequence_name=Y38C1AB.4"].extend(["IV","+",695,1,0])





# output the pormoter file
with open(promoters, 'w') as PROMOTERS:
	for ID,values in promoter_regions.items():
		chromosome,orient,start,promoter_start,distance=values
		if chromosome != "MtDNA":
			if orient == "+":
				start=start-1
				PROMOTERS.write("{chromosome}\tWormBase\tpromoter\t{promoter_start}\t{start}\t.\t{orient}\t.\t{ID}\n".format(**locals()))
			if orient == "-":
				start=start+1
				PROMOTERS.write("{chromosome}\tWormBase\tpromoter\t{start}\t{promoter_start}\t.\t{orient}\t.\t{ID}\n".format(**locals()))

# sort promoter file
os.system("cat {promoters}  | sort -k1,1 -k4,4n > tmp && mv tmp {promoters}".format(**locals()))




