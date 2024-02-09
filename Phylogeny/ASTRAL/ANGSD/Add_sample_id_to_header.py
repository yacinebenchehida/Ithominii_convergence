#!/usr/bin/python
####################
# Import libraries #
####################
from Bio import SeqIO
import sys
import os
import gzip 

############################################################
# load fasta sequence and the pattern to add to the header #
############################################################
sequence_fasta = sys.argv[1]
id =  str(sys.argv[2])

##############################
# Add pattern to the heaader #
##############################
with gzip.open(sequence_fasta, "rt") as handle:
	for record in SeqIO.parse(handle, "fasta"):
		print(">" + record.id + "_"  + str(id))
		print(record.seq)
