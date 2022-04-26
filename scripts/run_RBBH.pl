#!/usr/bin/perl

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 3 || $ARGV[0] eq "-h")
        {print "\nGiven two prot fasta file,find RBBH between them\n";
                print "Usage:\t $scriptname <fasta_file1> <fasta_file2> <species_name1> <species_name2>\n";
                print "Written by Guo qingxiang, guoqing\@webmail.hzau.edu.cn. \n";
                print "Distributed without any guarantees or restrictions\n";
                        print "\n"; exit;
                                }


open IN1, "$ARGV[0]";
open IN2, "$ARGV[1]";


$file1 = $ARGV[0];
$file2 = $ARGV[1];
$name1 = $ARGV[2];
$name2 = $ARGV[3]; 

system("blast_rbbh.py $file1 $file2 -o ${name1}_$name2.tab -t blastp -a prot -i 1 -c 1 --threads 32");
system("cat ${name1}_$name2.tab | cut -f 1,2 > tmp");
system("tail -n +2 tmp > ${name1}_$name2");
system("rm tmp");


