#!/usr/bin/perl

open IN1, "$ARGV[0]";

$number = $ARGV[1];

open OUT, ">out2";


while (<IN1>) {
	chomp;
	$_ =~ /(\S+)\s(\S+)\s(\S+)/;
	if ($3 gt $number) {
	print OUT "$_\n";

}


}
