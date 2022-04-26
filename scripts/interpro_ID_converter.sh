#!/bin/sh

grep "IPR" list  > out
split_lines_by_semicolon.pl out
grep "IPR" out2 > out3
remove_duplicate.pl out3
perl -p -i -e 's/\s+(\S+)$//g' duplicate_remove
mv duplicate_remove duplicate_remove1
remove_duplicate.pl duplicate_remove1
rm duplicate_remove1
mv duplicate_remove final_IPR
perl -p -i -e 's/^\s//g' final_IPR
remove_duplicate.pl final_IPR
rm final_IPR
mv duplicate_remove final_IPR
convert_interpro_2_name.pl final_IPR /home/train/public_database/Interpro/entry.list
rm out out2 out3
cat convert_name | cut -f 3 > domain_name


