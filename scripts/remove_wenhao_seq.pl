#!/usr/bin/perl

open IN, "$ARGV[0]";
open OUT, ">out";

while (<IN>) {
	chomp;
   if ($_ =~ /\?/) {

} else {
	print OUT "$_\n";

}



}
