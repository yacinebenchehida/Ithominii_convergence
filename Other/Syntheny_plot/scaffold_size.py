#!/usr/bin/python
from Bio import SeqIO
import sys
import os

input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), "fasta")
mid_point = round(float(sys.argv[2]))
gene = sys.argv[3]

for record in input_seq_iterator:
	if (gene == "Cortex1"):
		print(str(mid_point - 39654) + "\t" + str(mid_point + 39654))
	else:
		print(str(mid_point - len(record.seq)*10) + "\t" + str(mid_point + len(record.seq)*10))
