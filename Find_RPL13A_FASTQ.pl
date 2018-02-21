#!/usr/bin/perl
use strict;
#use warnings;
############################## 
my $ref="GATTGGACATTCGGAAGAGGGCCCGCCTTCCCTGGGGAATCTCTGCGCACGCGCAGAACGCTTCGACCAATGAAAACACAGGAAGCCGTCCGCGCAACCGCGTTGCGTCACTTCTGCCGCCCCTGTTTCAAGGGATAAGAAACCCTGCGACAAAACCTCCTCCTTTTCCAAGCGGCTGCCGAAGATGGCGGAGGTGCAGGTATGGGCTCCGCGCGGGCCGGGGCGGCAAGGGGCCGGGTGGGATCCAGG";
##############maybe need to make a foreach for all lanes
my $in= "/cs/cbio/oren/SequencingData/Hagar/ManarMiseq/B200_S1_L001_R1_001.fastq";
#mkdir $path unless(-d $path);
##########################################################################
open IN,"<$in";
my %hash;	
my $line=<IN>; chomp $line;
my $readcount=0;
my $allreadcount=0;
while (defined($line)){
	#if($line=~ m/^(\>.+)/){
	if($line=~ m/^\@(.+)/){
		#print $line."\n";
		my $head=$1;
		$line=<IN>; chomp $line;
		if (length($line)==249){
			$hash{$head}=$line;
			$readcount++;
			#print $line."\n";
		}
		$allreadcount++;
		$line=<IN>; chomp $line;
		$line=<IN>; chomp $line;
		$line=<IN>; chomp $line;
	}	
	
}
close(IN); 
print "allreadcount: ".$allreadcount."\n";
print "249readcount: ".$readcount."\n";
###########################################################
my $out= "/cs/cbio/oren/SequencingData/Hagar/ManarMiseq/Manar_matrix.txt";
#my $out1= "/cs/cbio/oren/SequencingData/Hagar/ManarMiseq/Manar_anti".$numrep.".txt";
#my $out2= "/cs/cbio/oren/SequencingData/Hagar/ManarMiseq/Manar_sense".$numrep.".txt";
#my $out3= "/cs/cbio/oren/SequencingData/Hagar/ManarMiseq/Manar_oops".$numrep.".txt";
#open OUT,">$out";
#open OUT1,">$out1";
#open OUT2,">$out2";
#open OUT3,">$out3";
my @ref=split("",$ref);
foreach my $key (keys(%hash)){
	my @b=split("",$hash{$key});
	my $score=0;
	my $g=0;
	my $c=0;
	my $a=0;
	my $t=0;
	my $gaind="";
	my $ctind="";
	for(my $x=0;$x<scalar(@ref);$x++){
		if ($b[$x] eq $ref[$x]){$score++;}
		if ($ref[$x] eq "G"){$g++;}
		if ($ref[$x] eq "C"){$c++;}
		if (($ref[$x] eq "G") && ($b[$x] eq "A")){$a++;$gaind=$gaind.$x.",";}
		if (($ref[$x] eq "C") && ($b[$x] eq "T")){$t++;$ctind=$ctind.$x.",";}
	}
	if ($gaind eq ""){$gaind="NA";}
	if ($ctind eq ""){$ctind="NA";}
	#print OUT $score."\t".$g."\t".$a."\t".$c."\t".$t."\t".$gaind."\t".$ctind."\n";
	#print $score."\t".$g."\t".$a."\t".$c."\t".$t."\t".$gaind."\t".$ctind."\n";
}
#close(OUT1); close(OUT2); close(OUT3);
#close(OUT); 
#print "number: ".$count."\n";
print "END\n";