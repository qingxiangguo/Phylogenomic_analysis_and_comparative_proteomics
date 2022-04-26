#!/usr/bin/perl

open IN, "$ARGV[0]";

$w = $ARGV[0];

open OUT, ">${w}_filter";


while (<IN>) {
	chomp;
	$_ =~ /(\S+)\t(\S+)/;

	our $part = $1;
	our $seq = $2;


	if ($part =~ /ENLEE|HESPP|KUIWA|MYAMP|MYCER|MYHON|MYPEN|MYPHY|MYTUR|MYWUL|POHYD|THKIT|SPZAH|TEBRY/) {
	
	print OUT "$part\n$seq\n";

} else {

}


}
