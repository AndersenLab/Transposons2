#!/usr/bin/env python

import re


abs_non_genes={}
abs_genes={}
ref_non_genes={}
ref_genes={}
ins_non_genes={}
ins_genes={}

non="dedup_nongene.txt"
with open(non, 'r') as IN:
	for line in IN:
		items = re.split("[\t]", line)
		method = items[1]
		event = items[3]
		geneID = items[7]
		if method == "absence":
			abs_non_genes[geneID]=line
		elif method == "insertion":
			ins_non_genes[geneID]=line
		elif method == "reference":
			ref_non_genes[geneID]=line
		else:
			print "error"


OUT=open("removed_genes.txt", 'w')

gen="dedup_genes.txt"
with open(gen, 'r') as IN:
		for line in IN:
			items = re.split("[\t]", line)
			method = items[1]
			event = items[3]
			geneID = items[7]
			if method == "absence":
				if geneID not in abs_non_genes.keys(): #only put in dictionary if the gene is not already represents as exon, into, UTR
					items[3] == "other" #reassign "gene" to "other"
					OUT.write(line)
			elif method == "reference":
				if geneID not in ref_non_genes.keys(): #only put in dictionary if the gene is not already represents as exon, into, UTR
					items[3] == "other" #reassign "gene" to "other"
					OUT.write(line)
			elif method == "insertion":
				if geneID not in ins_non_genes.keys(): #only put in dictionary if the gene is not already represents as exon, into, UTR
					items[3] == "other" #reassign "gene" to "other"
					OUT.write(line)
			else:
				print method
				print "error"

OUT.close()

