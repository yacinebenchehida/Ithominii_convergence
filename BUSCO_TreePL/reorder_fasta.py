#!/usr/bin/env python

import sys

#  reorder_fasta.py <file.fasta> <file.reference>

fasta= open(sys.argv[1])
ref= open(sys.argv[2])

seq_dict= {}
while True:
    line= fasta.readline()
    if line == '':
        break
    if line.strip().startswith('>'):
        seq_name= line.strip()[1:]
        seq_dict[seq_name]= []
    else:
        seq_dict[seq_name].append(line.strip())
fasta.close()
for seq_name in ref:
    seq_name= seq_name.strip()
    print('>' + seq_name)
    print('\n'.join(seq_dict[seq_name]))
ref.close()
sys.exit()
