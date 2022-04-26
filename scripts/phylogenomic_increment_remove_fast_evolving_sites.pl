#!/usr/bin/perl

$scriptname=$0; $scriptname =~ s/.+\///g;

if ($#ARGV != 2 || $ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "--help")
        {print "\nRead in a fasta sequence file and incrementally delete the fast evloving sites according to the result of Dist_Est\n";
        print "Usage:\t $scriptname <fasta_alignment> <rate_est.dat[output of Dist_Est]> <num_sites_per>\n";
        print "\n"; exit;
        }

open IN1, "$ARGV[0]";
open IN2, "$ARGV[1]";
open IN3, "$ARGV[2]";

our $num_sites = $ARGV[2];
our $length ;




while (<IN1>) {
	chomp;
	if ($_ =~ />/) {
} else {
	our $length = length ($_);   #get the length of alignment;
}
}
#---------prepare files----

 system ("cat $ARGV[1] | cut -d ' ' -f 4 > new1");
 system ("cat -n new1 > new2");
 system(q(perl -p -i -e "s/^\s+//g" new2));
 system ("sort -k2,2nr new2 > new3");
system ("flat_the_fasta_seq.pl $ARGV[0]");

#-------start the loop----

	our $count = 0;
	for ($i = $num_sites ; $i < $length ; $i += $num_sites)  {
	system ("head -n $i new3 | cat | cut -f 1 > new_$i");
	system("line_to_coma.sh new_$i");
	system ("phylogenomic_remove_specific_position_alignment.pl new_$i flated > out_${i}_removed.fasta");	
        system ("rm new_$i");
	$count ++ ;
}

	system ("rm new1 new2 new3 flated");


print "Your alignment length is <$length>\n\n\n";

print "You want to delete <$num_sites> fastest evolving sites per time\n\n\n";
print "Ok, I will produce <$count> deduced files for you\n\n\n";
