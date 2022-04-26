#!/bin/bash

for dir in ./*.fa 
do

NAME=`basename $dir \.fa`

remove_contaminant_by_ID.pl $dir ${NAME}\.out_best_delete

mv survive.fasta ${dir}_revised

done
