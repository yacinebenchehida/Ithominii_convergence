#!/bin/bash
#Author: Yacine Ben Chehida

SCAFFOLD=$1
WINDOW=$2
SCAF_SIZE=$3
SPECIES=$4
INPUT=$5

test -f $INPUT/$SPECIES/$SCAFFOLD/"$SCAFFOLD"_windows_position.txt  && rm $INPUT/$SPECIES/$SCAFFOLD/"$SCAFFOLD"_windows_position.txt

for ((i = 1, j = $WINDOW; i < $SCAF_SIZE && j < $SCAF_SIZE; i = j + 1, j=j+$WINDOW))
	do echo -e $i"\t"$j >> $INPUT/$SPECIES/$SCAFFOLD/"$SCAFFOLD"_windows_position.txt
done

j=$SCAF_SIZE
echo -e $i"\t"$j >> $INPUT/$SPECIES/$SCAFFOLD/"$SCAFFOLD"_windows_position.txt
