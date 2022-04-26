#!/bin/bash
#This shell script removes any sequences that are mostly gaps from a fasta alignment
#If more than $max_percent_gaps percent of the sequence is gaps, it is deleted
#In other words, $max_percent_gaps is the maximum allowed percentage of missing data / gaps in any one sequence.

#####CHANGE THIS VARIABLE AS NEEDED#####
max_percent_gaps=50
########################################

echo "Removing sequences that are more than $max_percent_gaps percent gaps..."

for FILENAME in *.fa
do
bashcode='sequ=\1;dashes=$(echo "$sequ" | tr -cd "-");echo $(( ${#dashes}*100/${#sequ} > '$max_percent_gaps' ))' #Determines if the alignment contains any sequences which are more than X% gaps/missing data
sed -r "/>/{N; /-/{h; s:.*\n(.*):$bashcode:e; /1/d; /0/x}}" $FILENAME > $FILENAME.long #Removes any sequences that are more than X% gaps/missing data and writes the output to $FILENAME.long
done

echo "Done!"

