#!/usr/bin/env python 
# old check script...ignore

import sys
import re
import os

IN_NEW = open("/lscr2/andersenlab/kml436/git_repos2/Transposons/files/new_WB_famlies.txt", "r")
family_renames={}
for line in IN_NEW:
	line = line.rstrip('\n')
	items= re.split("[\t]",line)
	old_name = items[0]
	new_name = items[1]
	family_renames[old_name] = new_name
IN_NEW.close()

in_file_info = sys.argv[1]
in_file = "{in_file_info}_telocate_nonredundant.bed".format(**locals())
print in_file

for key in family_renames.keys():
	value = family_renames[key]
	print key
	print value
	 # && mv test
	os.system("cat {in_file}  | sed 's/{key}/llllllll/g' > FF".format(**locals()))
