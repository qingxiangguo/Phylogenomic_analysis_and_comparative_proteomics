#!/usr/bin/perl


my @files = `ls $ARGV[0]`;


open CMDS, '>', "blast_cmds" or die $!;

foreach (@files) {
	chomp;
	next unless m/(\S+)\.fa_filter$/;
	print CMDS "blastp -query $ARGV[0]/$_ -db ./all -out $ARGV[0]/$1.out -outfmt 6 -evalue 1e-10 -max_hsps 1 -num_threads 2\n";
}

close CMDS;






