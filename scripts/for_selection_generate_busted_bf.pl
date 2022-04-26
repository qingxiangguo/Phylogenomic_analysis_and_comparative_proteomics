#!/usr/bin/perl

 open IN, "$ARGV[0]";

        while (<IN>) {
	chomp;
	my $name = $_;
	`generate_busted_hyphy.sh ../prot_guided_align/$name.fa.long ../prot_guided_align/$name.tre > ./$name.bf`;

}
