#!/bin/bash
MIN_TAXA=71
mkdir -p rejected_few_taxa_2
for FILENAME in *.trim
do 
grep ">" $FILENAME > $FILENAME\.tmp
cat $FILENAME\.tmp | cut -f 1 -d "|" > $FILENAME\.tmp2
uniq $FILENAME\.tmp2 > $FILENAME\.tmp3
taxon_count=`grep -v 0 $FILENAME\.tmp3 | wc -l`

if [ "$taxon_count" -lt "$MIN_TAXA" ] ; then
echo $FILENAME
mv $FILENAME ./rejected_few_taxa_2/
fi
done

rm -rf *.tmp *.tmp2 *.tmp3
