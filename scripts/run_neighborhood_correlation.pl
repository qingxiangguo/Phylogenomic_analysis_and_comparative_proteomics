#!/usr/bin/perl

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 4 || $ARGV[0] eq "-h")
        {print "\nGiven two prot fasta file, using NC method to find homology between them\n";
                print "Usage:\t $scriptname <fasta_file1> <fasta_file2> <species_name1> <species_name2> <NC_threshold>\n";
                print "Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n";
                print "Distributed without any guarantees or restrictions\n";
                        print "\n"; exit;
                                }


open IN1, "$ARGV[0]";
open IN2, "$ARGV[1]";
 $nc_number = $ARGV[4];


$file1 = $ARGV[0];
$file2 = $ARGV[1];
$name1 = $ARGV[2];
$name2 = $ARGV[3]; 

system("cat $file1 $file2 > all.fasta");
system("makeblastdb -in all.fasta -dbtype prot -out all");
$seq_number = `grep ">" all.fasta -c`;
system("fastalength all.fasta > stat");
system("cat stat | cut -d ' ' -f 1 > w");
$residue_number = `awk '{ sum += \$1 } END { print sum }' w`;

$e = $seq_number * 10;
$searchsp = $residue_number * $residue_number;

system("blastp -db all -query all.fasta -out tab -evalue $e -num_threads 48 -outfmt 6 -searchsp $searchsp");
system("cat tab | cut -f 1,2,12 > raw");
system(q(perl -p -e "s/\t/ /g" raw > dat));
system("NC_standalone -f dat -o out --num_residues $residue_number");
system("NC_thresh.pl out $nc_number");
system("proteome_NC_result_process.pl out2");
system("mv NC_result ${name1}_$name2");
system("rm w stat raw out2");
