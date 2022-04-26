#!/usr/bin/perl

 open IN1, "$ARGV[0]";
open  OUT, ">partition_data.txt";

   while (<IN1>) {

	$_ =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/;
	
	print OUT "AUTO, $1 = $2-$3\n";


}
