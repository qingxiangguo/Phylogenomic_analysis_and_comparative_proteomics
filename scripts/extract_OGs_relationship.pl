#!/usr/bin/perl

open IN, "$ARGV[0]";
open OUT, ">rela";
while (<IN>) {

	chomp;
	
	if ($_ =~ /^(\d+?)\|\S+\s+(\d+?)\|/) {
	print OUT "$1\t$2\n";


}

}
