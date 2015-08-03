#!/usr/bin/env python 
# this script takes the chromosome regions to be available for the RSVSim analysis and outputs them in the appropriate format 
# to be used for the RSVSim simulations
# USE: prep_RSVSim_list.py

import re
import os
from  collections import defaultdict

positions_to_keep = {}
positions_to_keep = defaultdict(list)



in_file = "/lscr2/andersenlab/kml436/git_repos2/Transposons/files/RSVSim.bed"
IN_FILE= open(in_file, "r")

for line in IN_FILE:
	items = re.split("[\t]", line)
	chromosome = items[0]
	chromosome = str(chromosome)
	start_interval = int(items[1])
	end_interval = int(items[2])
	interval = "{start_interval}:{end_interval}".format(**locals())
	print interval
	positions_to_keep[chromosome].append(interval)
	#		print "removed"
	#		print i
out ="RSVSim_chromosomes.txt"
OUT_FILE = open(out, "w")

#average_FDR[key] = map(s, average_FDR[key]) 
for i in positions_to_keep.keys():
	value = positions_to_keep[i]
	print value
	OUT_FILE.write("{i}\tc(".format(**locals()))
	for v in value:
		print v
		OUT_FILE.write("{v},".format(**locals()))
	OUT_FILE.write(')\n')

IN_FILE.close()
OUT_FILE.close()
os.system("cat RSVSim_chromosomes.txt| sed 's/,)/)/g' > temp && mv temp RSVSim_chromosomes.txt")
	



