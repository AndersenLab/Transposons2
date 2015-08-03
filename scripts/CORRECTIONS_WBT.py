#!/usr/bin/env python 
# old check script...ignore

import sys
import re
import os

IN_NEW = open("/lscr2/andersenlab/kml436/git_repos2/Transposons/files/new_WB_familes.txt", "r")
names={}
for line in IN_NEW:
	line = line.rstrip('\n')
	items= re.split("[\t]",line)
	old_name = items[0]
	new_name = items[1]
	names[old_name] = 0
IN_NEW.close()

wbts="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/WBTs.txt"
WBTS= open(wbts, "r")

for line in WBTS:
	line = line.rstrip('\n')
	org_name = line
	if org_name not in names.keys():
		print org_name
	
WBTS.close()




