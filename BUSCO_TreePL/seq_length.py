#!/usr/bin/python
from Bio import SeqIO
import sys
import os

#usage: python -W ignore *fasta

input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), "fasta")
indels = 0 
tot_size = 0

for record in input_seq_iterator:
	sequences = str(record.seq)
	indels = indels + sequences.count('-')
	tot_size = tot_size + len(record.seq)
	pct_indels = (indels/tot_size)*100

print("{:.0f}".format(len(record.seq)))
