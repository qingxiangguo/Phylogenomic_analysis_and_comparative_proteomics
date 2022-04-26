#!/usr/bin/perl
use strict;
my $usage = <<USAGE;
Usage:
	perl $0 SingleCopyOrthologGroups 24[this is CPU threads]

USAGE
die $usage if @ARGV != 2;

my @files = `ls $ARGV[0]`;

open CMDS, '>', "fasttree_cmds" or die $!;

foreach (@files) {
	chomp;
	next unless m/(\S+)\.fa\.long$/;
	print CMDS "FastTreeMP -slow -gamma $ARGV[0]/$_ > $ARGV[0]/$1.tre\n";
}

close CMDS;

my $parafly_result =`ParaFly -c fasttree_cmds -CPU $ARGV[1] -shuffle -v`;

#unlink "clustalo_cmds", "clustalo_cmds.completed" unless -e "FailedCommands";

