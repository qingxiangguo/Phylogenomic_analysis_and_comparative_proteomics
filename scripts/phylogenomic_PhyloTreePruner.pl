#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
	perl $0 SingleCopyOrthologGroups 24[this is CPU threads]

USAGE
die $usage if @ARGV != 2;

my @files = `ls $ARGV[0]`;

open CMDS, '>', "PhyloTreePruner_cmds" or die $!;

foreach (@files) {
	chomp;
	next unless m/(\S+)\.fa\.long$/;
	print CMDS "java -cp /usr/local/bin/ PhyloTreePruner $ARGV[0]/$1.tre 50 $ARGV[0]/$_ 0.5 r\n";
}

close CMDS;

my $parafly_result =`ParaFly -c PhyloTreePruner_cmds -CPU $ARGV[1] -shuffle -v`;

#unlink "clustalo_cmds", "clustalo_cmds.completed" unless -e "FailedCommands";

