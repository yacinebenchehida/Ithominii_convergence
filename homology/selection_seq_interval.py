#!/usr/bin/python
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os

#usage: python -W ignore *fasta contigs start end 

input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), "fasta")
start = int(sys.argv[3]) - 1 
end =  int(sys.argv[4]) 
sp = sys.argv[5]

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp_sequence = ''.join(complement_dict[base] for base in reversed(dna_sequence))
    return reverse_comp_sequence


for record in input_seq_iterator:
        if(record.id == sys.argv[2]):
                if(sp == "Hypothyris_anastasia"):
                        print(">"+ str(record.id))
                        my_seq=record.seq[start:end]
                        print(reverse_complement(my_seq.upper()))
                else:
                        print(">"+ str(record.id))
                        print(record.seq[start:end])
