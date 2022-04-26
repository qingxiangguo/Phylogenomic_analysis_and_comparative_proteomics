#!/usr/bin/perl


$name1 = $ARGV[1];
$name2 = $ARGV[2];

open IN, "$ARGV[0]";
open OUT, ">${name1}_${name2}";

while (<IN>) {
	chomp;
	if ($_ =~ /$name1/) {
	if ($_ =~ /$name2/) {
	print OUT "$_\n";
}
}
}
