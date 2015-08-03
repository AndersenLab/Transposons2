#!/usr/bin/bash
# runs trimal
# USE:run_trimal.sh
# example:/lscr2/andersenlab/kml436/git_repos2/Transposons/files/alignments$ bash ../../scripts/run_trimal.sh
# ignore this script
#filename="$1"

while read -r line
do
    echo "Name read from file - $line"
    arr=($line)
    file_name=${arr[0]}
    parameters=${arr[1]}
    echo $file_name
    echo $parameters
    name=$(basename "$file_name" .aln)
    #adjust positions
    /lscr2/andersenlab/kml436/git_repos2/trimal/source/trimal -in $file_name -out EDITED_${file_name} -selectcols { $parameters } 
    /lscr2/andersenlab/kml436/git_repos2/Sequence-manipulation/Consensus.pl -in EDITED_${file_name}  -out CONSENSUS_${name}
    cat CONSENSUS_${name} ${name} > check_consensus_${name}
    /opt/mafft-6.240/bin/mafft check_consensus_${name} > A_check_consensus_${name}.aln
    cp CONSENSUS_${name} ../final_consensus_fastas

done < "trim_info.txt"

#simplify names in final consensus directory
cd ../final_consensus_fastas
rename 's/^CONSENSUS_//' *.fasta
rename 's/_fastas//' *.fasta