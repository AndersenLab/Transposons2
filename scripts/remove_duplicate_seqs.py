#!/usr/bin/env python 
# this script takes a fasta file as input and removes seuquences that are exact duplciates of one another
# for duplicates, one representative seqeunce will be saved to the ouput file and the names of
# all sequential duplicates skipped will be outut to a log file
# outout file with no duplicates ends in "_no_dups_fasta"
# USE: remove_duplicate_seqs.py <fasta_file>
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
import sys
import os

all_seqs = sys.argv[1]
file_name=os.path.splitext(all_seqs)[0]
out=file_name+ str("_no_dups.fasta")
log_out ="duplicate_log_"+str(file_name)+str(".log")
LOG_OUT = open(log_out, "w")
def remove_dup_seqs(records):
    """"SeqRecord iterator to removing duplicate sequences."""
    checksums = set()
    for record in records:
        checksum = seguid(record.seq)
        if checksum in checksums:
            LOG_OUT.write("Ignoring %s \n" % record.id)
            continue
        checksums.add(checksum)
        yield record

records = remove_dup_seqs(SeqIO.parse(all_seqs, "fasta"))
count = SeqIO.write(records, out, "fasta")
print "Saved %i records" % count
