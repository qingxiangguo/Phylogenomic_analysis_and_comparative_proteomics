#!/usr/bin/perl
use IO::File;
open IN1, "$ARGV[0]";


system('mkdir OGs');

while (<IN1>) {
	chomp;
$_ =~ /^(\d+)<=>(\S+)/;

    $number = $1;
my $fh=IO::File->new("> ./OGs/$number.txt");

	print $fh "$2\n";

}


