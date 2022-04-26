#!/usr/bin/perl

open IN1, "$ARGV[0]";  #group.txt

open OUT1, ">header";


while (<IN1>) {
        chomp;
	$_ =~ /^(\S+):\s(.*)/;
	$ocg =$1;
	$left = $2;
	my @fields = split /\s+/, $2 ;

	foreach $field (@fields) {
	print OUT1 "$ocg\t$field\n";



}


}



