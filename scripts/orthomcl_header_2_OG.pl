#!/usr/bin/perl

open IN1, "$ARGV[0]";  #seq
open IN2, "$ARGV[1]";  #header

$name = $ARGV[0];


open OUT1, ">${name}_OG";

chomp(@lines = <IN1>);


while (<IN2>) {
  	   chomp;
	$_ =~ /(\S+)\t(\S+)/;
	$hash{$2} = $1 ;

}



foreach $lines (@lines) {
        if (exists $hash{"$lines"}) {
        print OUT1 "$hash{$lines}\n";
} else {
        print OUT1 "$lines\n";

}
}



	




