#!/bin/bash
# this script reformats the repbase classifications to remove spaces forward slashes from the transposon classification  names
# USE: reformat_repbase_headers.sh


sed 's/DNA transposon/DNA_transposon
sed 's/EnSpm\/CACTA/EnSpm_CACTA
sed 's/LTR Retrotransposon/LTR_Retrotransposon
sed 's/Mariner\/Tc1/Mariner_Tc1
sed 's/SINE2\/tRNA/SINE2\/tRNA
sed 's/Transposable Element/