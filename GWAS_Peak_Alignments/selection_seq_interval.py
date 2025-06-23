#!/usr/bin/python
from Bio import SeqIO
from Bio.Seq import Seq
import sys

# Usage: python script.py input.fasta contig_name start end species_name

input_fasta = sys.argv[1]
contig_name = sys.argv[2]
start = int(sys.argv[3]) - 1  # Convert to 0-based index
end = int(sys.argv[4])
species_name = sys.argv[5]

sign = end - start

if(sign < 0):
        tmp = start
        start = end
        end = tmp

# Function to compute reverse complement using Biopython
def reverse_complement(dna_sequence):
	return str(Seq(dna_sequence).reverse_complement())

# Read the input FASTA file
input_seq_iterator = SeqIO.parse(open(input_fasta, "rU"), "fasta")

# Iterate through sequences
for record in input_seq_iterator:
	if record.id == contig_name:
		sequence = str(record.seq)
		selected_sequence = sequence[start:end]

		# Check if reverse complement is needed based on species name
		if species_name in ("Hypothyris_anastasia", "Melinaea_mothone"):
			print(f">{record.id}")
			print(reverse_complement(selected_sequence.upper()))
		else:
			print(f">{record.id}")
			print(selected_sequence)
