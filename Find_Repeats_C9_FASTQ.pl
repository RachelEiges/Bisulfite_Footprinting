#!/usr/bin/perl
use strict;
#use warnings;
#####################CAT
#>M01015:70:000000000-B45YW:1:1101:16264:1748 2:N:0:ATCTCAGG+TCTGCATA
#TACCTGATAAAGATTAACCAGAAGAAAACAAGGAGGGAAACAACCGCAGCCTGTAGCAAGCTCTGGAACTCAGGAGTCGCGCGCTAGGGGCCGTGGCCGTGGCGTGGTCGGGGCGGGCCCCAAACCAAACCCCAGCCCAACCTTCGGTTTCGGTACCTTCCACCGCGGCGTGCGAGGCCCCAGTGGGTGGGTGTGTGGGTGCGGGGGGTGGGCGTCTCGGCGGGGGGGGTGTTTGTGTGGTTGTTGGGGG
#>M01015:70:000000000-B45YW:1:1101:16934:1748 2:N:0:ATCTCAGG+TCTGCATA
#GACCTGATAAAGATTAACCAGAAGAAAACAAGGAGGGAAACAACCGCAGCCTGTAGCAAGCTCTGGAACTCAGGAGTCGCGCGCTCGTGGCCGGGGCCGGGGCGTCGTCGAAACCACCCCAGGGACAGACCCGAGGCGGGGCTGCGGTTGTGGGGCCTTCGCCCTCCGCGGGGGGGGCGCGGGGGGGGGGGGGGGTGTGGGGGGGGGGGGGGGCGCTCTGGGTGGGGGGGGTGTGGGGGTGGGGGGGGGG
#>M01015:70:000000000-B45YW:1:1101:15399:1760 2:N:0:ATCTCAGG+TCTGCATA
#TACCTGATAAAGATTAACCAGAAGAAAACAAGGAGGGAAACAACCGCAGCCTATAGCAAGCTCTGGAACTCAGGAGTCGCGCACTAGGGGCCGTGGCCAAAACTTCATCAAAACAAACCCAAAAACCAACCCAGGGCGGGACTGCGGTTGCGGTGCCTGCCCCCACGTCGGCGCATCCGCAGGCGGTGGGGTGGGTGTGGGTGTCGGGTGGGGCTCTCGTGGGGGGGTGGTTTTGTGGTGGGTGGTTGGG
#>M01015:70:000000000-B45YW:1:1101:13357:1763 2:N:0:ATCTCAGG+ACTACATA
#TACCTGATAAAGATTAACCTGAAGAAAACAAGGAGGGAAACAACCGCAGCCTTTAGCAAGCTCTGGAACTCAGGAGTCGCGCGCTAGGGGCCGTGGCCGGGGCTTGGTCGGAGCGCGCCCAACAACCACCCCCACCCAAGCCTGCGGTTGCGGTGCCTGCGCCCGCGGCGGCGTCGTCGCAGGCGGTGGCGGGGTGGCTCGTTGGGGGGGGGCGTCCTTCGGGTGGGGTGGTTTGCGGTCGGCTGCGGGG
#GGGGCCGGGGCC     
############################## 
##############maybe need to make a foreach for all lanes
my $in= "/cs/cbio/oren/SequencingData/Hagar/C9/Hagar3_ALS3_2_L001_R2_001.gz.fa";
#mkdir $path unless(-d $path);
##########################################################################
open IN,"<$in";
my %hash;	
my $line=<IN>; chomp $line;
my $readcount=0;
my $allreadcount=0;
while (defined($line)){
	if($line=~ m/^\>(.+)/){
	#if($line=~ m/^\@(.+)/){
		#print $line."\n";
		my $head=$1;
		$line=<IN>; chomp $line;
		if (length($line)==250){
			$hash{$head}=$line;
			$readcount++;
			#print $line."\n";
		}
		$allreadcount++;
		$line=<IN>; chomp $line;
		#$line=<IN>; chomp $line;
		#$line=<IN>; chomp $line;
	}	
	
}
close(IN); 
print "allreadcount: ".$allreadcount."\n";
print "250readcount: ".$readcount."\n";
my @anti=("AGGGCC","GAGGCC","GGAGCC","GGGACC","AAGGCC","AAAGCC","AAAACC","GAAGCC","GGAACC","GAAACC","AAGACC","AGAACC","AGGACC","AGAGCC","GAGACC");
my @sense=("GGGGTT","GGGGCT","GGGGTC");
foreach my $numrep (2..2){
	print "number of repeats is: ".$numrep."\n";
	#my $out= "/cs/cbio/oren/SequencingData/Hagar/Hagar4_SZ13.fa";
	my $out1= "/cs/cbio/oren/SequencingData/Hagar/C9/Hagar3_ALS3_anti".$numrep.".txt";
	my $out2= "/cs/cbio/oren/SequencingData/Hagar/C9/Hagar3_ALS3_sense".$numrep.".txt";
	my $out3= "/cs/cbio/oren/SequencingData/Hagar/C9/Hagar3_ALS3_oops".$numrep.".txt";
	#open OUT,">$out";
	open OUT1,">$out1";
	open OUT2,">$out2";
	open OUT3,">$out3";
	my $oopscount=0;
	my $countA=0;
	my $countS=0;
	my $countO=0;
	my @count;
	$count[0]=0;
	foreach my $n (1..($numrep-1)){
		$count[$n]=$count[$n-1]+6;
		#print $count[$n]."\n";
	}
	foreach my $key (keys(%hash)){
		print OUT ">".$key."\n".$hash{$key}."\n";
		#print $hash{$key}."\n";
		my $a=substr($hash{$key}, 86, ($numrep*6));
		#my $a=substr($hash{$key}, 83, 60);
		my @b=split("",$a);
		my @c;
		#my @count=(0,3,6,9,12,15,18,21,24,27,30,33);
		#my @count=(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87);
		#my @count=(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57);
		my $anticount=0;
		my $sensecount=0;
		my $flagAnti=0; my $flagSense=0;
		foreach my $x (0..($numrep-1)){
			$c[$x]=$b[$count[$x]].$b[$count[$x]+1].$b[$count[$x]+2].$b[$count[$x]+3].$b[$count[$x]+4].$b[$count[$x]+5];
			#print $c[$x]."\t";
			if ($c[$x] eq "GGGGCC"){$anticount++;$sensecount++;}
			for (my $n=0;$n<scalar(@anti);$n++){
				if ($c[$x] eq $anti[$n]){$anticount++; $flagAnti=1;}
			}
			for (my $n=0;$n<scalar(@sense);$n++){
				if ($c[$x] eq $sense[$n]){$sensecount++; $flagSense=1;}
			}
		}
		if (($flagAnti==1) && ($flagSense==0)){$countA++; print OUT1 ">".$key."\n".$hash{$key}."\n";}
		elsif (($flagAnti==0) && ($flagSense==1)){$countS++; print OUT2 ">".$key."\n".$hash{$key}."\n";}
		else{print OUT3 ">".$key."\n".$hash{$key}."\n"; $countO++;}	
	}
	close(OUT1); close(OUT2); close(OUT3);
	print "number A: ".$countA."\n";
	print "number S: ".$countS."\n";
	print "number O: ".$countO."\n";
	#print "oopscount: ".$oopscount."\n"."END\n";
}
