#!/bin/bash
# this script transposes the rows and columns in a file
# USE: transpose_matrix.sh <input_file> <name_output_file>

input=${1}
output=${2}

awk '{for (i=1; i<=NF; i++) a[i]=a[i](NR!=1?FS:'\t')$i} END {for (i=1; i in a; i++) print a[i]}' $input  |sed 's/ /\t/g' > $output
