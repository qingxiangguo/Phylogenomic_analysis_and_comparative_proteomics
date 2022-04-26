#!/usr/bin/perl
	
$file = $ARGV[0];
open IN1, "$ARGV[0]";

$wd = `readlink -f $file`;

 chomp($wd);

	my ($DIR) = $wd =~ /\/(\w+)\.txt$/;

	print $DIR;

open OUT, ">${DIR}.fasta";

while (<IN1>) {
	chomp;
    $_ =~ /^(\S+\w\|)(\S+$)/;

	print OUT ">$1\n$2\n";

}


