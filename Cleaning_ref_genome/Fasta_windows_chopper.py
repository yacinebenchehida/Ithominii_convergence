#!/usr/bin/env python
import itertools
import sys

with open(sys.argv[4]) as f:
	for line1 in itertools.zip_longest(*[f]*2):
		line1 = [x.strip() for x in line1]
		if(line1[0].__contains__(sys.argv[3])):
			print(line1[0])
			print(line1[1][int(sys.argv[1])-1:int(sys.argv[2])])
