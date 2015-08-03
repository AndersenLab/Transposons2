#!/usr/bin/env python 
# this script goes through the telocate output and only saves calls in which the reference te and called te match and makes sure the 
# same reference te is not counted for more than once (will save information of call with highest read support)
# USE:process_telocate_output.py <closest_prepped.txt>

import sys
import re

te_calls=sys.argv[1]
TE_CALLS = open(te_calls, "r")

out = "processed_calls.txt"
OUT = open(out, "w")

TE_dict={}
for line in TE_CALLS:
	items = re.split("[\t]",line)
	chromosome = items[0]
	start = items[1]
	end = items[2]
	reference_TE = items[3]
	ID_reference = "{chromosome}_{start}_{end}_{reference_TE}".format(**locals())
	orient = items[5]
	call = items[9]
	read_support = items[10]

	match=re.search("(.*)_(\w+-)?reference", call)
	called_TE = match.group(1)

	if reference_TE == called_TE:
		if ID_reference not in TE_dict.keys():
			#want coordinates of the reference but all the other info from the call
			TE_dict[ID_reference]="{chromosome}\t{start}\t{end}\t".format(**locals()) + "\t".join(items[9:11]) + "\t{orient}\n".format(**locals())
		else:
			previous = re.split("[\t]", TE_dict[ID_reference])
			prev_read_support = previous[4]
			if read_support > prev_read_support: # replace if the next read support is higher
				TE_dict[ID_reference]="{chromosome}\t{start}\t{end}\t".format(**locals()) + "\t".join(items[9:11]) + "\t{orient}\n".format(**locals())
for key in TE_dict.keys():
	value = TE_dict[key]
	OUT.write(value)
