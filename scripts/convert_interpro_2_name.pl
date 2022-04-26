#!/usr/bin/perl

open IN1, "$ARGV[0]"; #list
open IN2, "$ARGV[1]"; #entry

open OUT1, ">convert_name";


chomp(@lines = <IN1>);


while (<IN2>) {
        chomp;
                $_ =~ /^(\S+)\t(.*)$/;

                        $hash_1{$1} = $_;

                        }



foreach $lines (@lines) {
        if (exists $hash_1{"$lines"}) {
        print OUT1 "$hash_1{$lines}\n";
} else {
	print OUT1 "0\n";

}
}
