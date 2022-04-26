#!/usr/bin/perl

open IN1, "$ARGV[0]";  #ab
open IN2, "$ARGV[1]";  #ac
open IN3, "$ARGV[2]";  #bc

open OUT1, ">RBBH_intersection";

while (<IN1>) {
	chomp;
	$_ =~ /(\S+)\t(\S+)/;
	$hash_a{$2} = $1;

}


while (<IN2>) {
	chomp;
	$_ =~ /(\S+)\t(\S+)/;
        $hash_b{$2} = $1;


}


while (<IN3>) {
	chomp;
	$_ =~ /(\S+)\t(\S+)/;

	if ($hash_a{$1} eq $hash_b{$2} && exists $hash_a{$1} && exists $hash_b{$2}) {

	print OUT1 "$_\t$hash_a{$1}\n";

}



}
