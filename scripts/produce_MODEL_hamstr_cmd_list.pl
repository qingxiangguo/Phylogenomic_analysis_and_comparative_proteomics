#!/usr/bin/perl

$file = $ARGV[0];


      $wd = `readlink -f $file`;
	chomp($wd);

	my ($DIR) = $wd =~ /\/(\w+)$/;	
	
open OUT, ">${DIR}_cmd.sh";

             $name1    = `ls $file/*fasta`;  
	   my  ($name2) = $name1 =~ /\S+\/(\w+)\.fasta/;  
          
print OUT "#!/bin/bash\n";

print OUT "cd $wd/model\n";

print OUT "nohup hamstr -protein -sequence_file=../$name2.fasta -taxon=$name2 -hmmset=modelorganisms_hmmer3 -refspec=DROME -relaxed -central -force -cpu 6 -eval_blast=0.01 -eval_hmm=0.01 &\n";


system ("chmod 755 ./${DIR}_cmd.sh");
