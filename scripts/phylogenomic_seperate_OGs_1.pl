#!/usr/bin/perl

open IN1, "$ARGV[0]";
open OUT, ">out";


while (<IN1>) {
	chomp;
	$_ =~ /^(\d+)\|/;
	$OG = $1;

	$hash{$_} = $OG ;

}


foreach $key(sort keys %hash){
if(exists $hash1{$hash{$key}}){
$hash1{$hash{$key}}.=",$key";}
else{
$hash1{$hash{$key}}=$key;
}
}


foreach $key(sort keys %hash1){
print OUT  "$key<=>$hash1{$key}\n";
}
