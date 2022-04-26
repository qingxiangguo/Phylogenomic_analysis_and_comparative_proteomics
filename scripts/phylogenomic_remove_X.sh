#!/bin/bash

for FILENAME in *.fasta
do
sed 's/^[^>]\{,19\}X//' $FILENAME > $FILENAME.trim
done
rename .fasta.trim .fasta *.fasta.trim

for FILENAME in *.fasta
do
sed '/>/! s/X.\{,19\}$//' $FILENAME > $FILENAME.trim
done
rename .fasta.trim .fasta *.fasta.trim

1
