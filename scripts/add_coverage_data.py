#!/usr/bin/env python 
# this script combines the information on total te counts and average coverage for a sample into one file
# USE: add_coverage_data.py <Full_Results> <coverage/bam_stats_file>
import sys
import re
from collections import defaultdict

results_file=sys.argv[1]
RESULTS_FILE=open(results_file,"r")
coverage_file=sys.argv[2]
COVERAGE_FILE=open(coverage_file, "r")
SAMPLE_TES_AND_COVERAGE=open("coverage_and_te_counts.txt", "w")

first_line=True
coverage={}
totals={}
totals=defaultdict(list)
for line in RESULTS_FILE:
	line=line.rstrip('\n')
	items=re.split("[\t]",line)
	if first_line:
		total_absent=items.index('absent_TRANS_total') #get index number of column name
		total_new=items.index('new_TRANS_total')
		total_reference=items.index('reference_TRANS_total')
	else:
		sample=items[0]
		totals[sample].extend((items[total_absent],items[total_new],items[total_reference])) # use extend to append multiple items to a list
	first_line=False
for line in COVERAGE_FILE:
	line=line.rstrip('\n')
	items=re.split("[\t]",line)
	bam=items[0]
	field=items[2]
	statistic=items[3]
	stat_value=items[4]
	if field =="BAM Statistics - Merged" and statistic=="Depth of Coverage (genome)":
		stat_value=round(float(stat_value),2)
		coverage[bam]=stat_value
		print stat_value
SAMPLE_TES_AND_COVERAGE.write("sample\tabsence\tinsertion\treference\tcoverage\n")
#remove the 2 unnecessary QX strains
#del totals['QX2265']
#del totals['QX2266']
for i in sorted(totals.keys()):
	SAMPLE_TES_AND_COVERAGE.write(i)
	for value in totals[i]:
		SAMPLE_TES_AND_COVERAGE.write("\t")
		SAMPLE_TES_AND_COVERAGE.write(value) #the sample name
	SAMPLE_TES_AND_COVERAGE.write("\t")
	SAMPLE_TES_AND_COVERAGE.write(str(coverage[i])) #the coverage level
	SAMPLE_TES_AND_COVERAGE.write("\n")

