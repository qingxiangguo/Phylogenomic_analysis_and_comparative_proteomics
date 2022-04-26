#!/usr/bin/perl

 open IN1, "$ARGV[0]";
open  OUT, ">charset";

   while (<IN1>) {

	$_ =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/;
	
	print OUT "charset $1 = $2 - $3 ; \n";


}
