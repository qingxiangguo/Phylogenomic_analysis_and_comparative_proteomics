#!/usr/bin/perl

open IN, "$ARGV[0]";

chomp(my @lines = <IN>);

foreach $line (@lines) {

	system ("rm ./${line}.fasta");


}
