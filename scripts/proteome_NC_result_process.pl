#!/usr/bin/perl

open IN1, "$ARGV[0]"; #union
open IN2, "$ARGV[0]"; #union

open OUT1, ">NC_result";

our $name;

while (<IN1>) {
	chomp;
	$_ =~ /^(\d+?)_(\S+)\s/;
	our	$name = $2;
	last;
}


close IN1;

while (<IN2>) {
        chomp;
	
	$_ =~ /(\S+)\s(\S+)\s(\S+)/;

	$name1 = $1;
	$name2 = $2;
	$name3 = $3;

	if ($name1 =~ /$name/ && $name2 !~ /$name/ ) {

	print OUT1 "$_\n";

}




                        }

