#!/bin/bash

MIN_SEQUENCE_LENGTH=100

for FILENAME in *.fasta
do
grep -B 1 "[^>].\{$MIN_SEQUENCE_LENGTH,\}" $FILENAME > $FILENAME.out
sed -i 's/--//g' $FILENAME.out
sed -i '/^$/d' $FILENAME.out
rm -rf $FILENAME
mv $FILENAME.out $FILENAME
done

