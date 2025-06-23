#!/usr/bin/python
from Bio import SeqIO
import sys
import os

#usage: python -W ignore *fasta contigs start end 

input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), "fasta")
start = int(sys.argv[3]) - 1 
end =  int(sys.argv[4])

for record in input_seq_iterator:
	if(record.id == sys.argv[2]):
		print(">"+ str(record.id))
		print(record.seq[start:end])
