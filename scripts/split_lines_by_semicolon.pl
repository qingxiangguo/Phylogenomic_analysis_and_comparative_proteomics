#!/usr/bin/perl

	open IN1, "$ARGV[0]";
	open OUT, ">out2";

while (<IN1>) {
	   	chomp;
             @array = split /;/, "$_";

	foreach $array (@array) {

	print OUT "$array\n";


}



}
