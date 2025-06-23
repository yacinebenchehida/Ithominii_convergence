#!/usr/bin/env python
from Bio import SeqIO
import sys
import os


input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), "fasta")
for record in input_seq_iterator:
	fasta_to_write = sys.argv[2] + "/" + record.id + "/" + record.id + ".fasta"
	output_handle = open(fasta_to_write, "w")
	SeqIO.write(record, output_handle, "fasta-2line")
	output_handle.close()
