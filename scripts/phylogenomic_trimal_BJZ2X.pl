#!/usr/bin/perl

	open IN1, "$ARGV[0]";
	$name =$ARGV[0];
	open OUT, ">$name.BJZ2X";
while (<IN1>) {
	chomp;
	if ($_ =~ /^>/) {


} else  {
	s/B|J|Z/X/g;
	

}

	print OUT "$_\n";

}
