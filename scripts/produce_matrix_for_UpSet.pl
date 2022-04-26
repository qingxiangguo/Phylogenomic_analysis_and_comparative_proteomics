#!/usr/bin/perl


  open IN1, "$ARGV[0]"; #gene
  open IN2, "$ARGV[1]"; #name

$name =  $ARGV[1];

  open OUT, ">${name}_matrix";

 while (<IN2>) {
	chomp;
	$hash{$_} = 1;

}


while (<IN1>) {
	chomp;
	if (exists $hash{$_}) {

	print OUT "1\n";
} else {
	print OUT "0\n";

}


}

