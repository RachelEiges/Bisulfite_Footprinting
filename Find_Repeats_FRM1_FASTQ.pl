#!/usr/bin/perl
use strict;
#use warnings;
#####################CAT
#>M01015:70:000000000-B45YW:1:1101:18781:1795 2:N:0:ACTCGCTA+AAGGAGTA
#CACCAGCTCCTCCATCTTCTCTTCAGCCCTGCTAGCGCCGGGAGCCCGCCCCCGAGAGGTGGGCTGCGGGTGTTTGAGGTCCAGCCGCCGCCGCCCGCGCTTGCCGTCGGCCCGTCGCCCGCTCAGAGGCGGCCCTCCACCGGAAGTGACACCGAACCGGGGCTGAGCGCCTCACTGCGGCCCCACCCCCGGCCCCCTTCGGGGGTGAACCCTGAAACCACGTCCCGCGCTCCACGCGGTTCCCTGTC
#####################
#@M01015:81:000000000-BCPFG:1:1101:15490:1688 2:N:0:CTCTCTAC+GAGCCTTA
#CACCAGCTCCTCCATCTTCTCTTCAGCCCTGCTAGCGCCGGGAGCCCGCCCCCGTGAGGCCGCCGCGCTGCCGCACGCCCCCTGGCAGCGGCGCCTCCGTCACCGCCGCTGCCCGCGCTCGCCGTCGGCCCGCCGCCCGCTCAGAGGCGGCCCTCCCCCGGCAGTGCAACCGAAACGGCGCTGAGCGCCTGTCTTCGGCCGACCCCCCGGCCTGCTGCGGGTGTTATCACTGACACCACTTCACCTGAT
#+
#11>11BAC11B>B1BAFGGGFGHB31ABAGA0G111000/////AAF/AEEEE/////0>F//>/>/>E/<F/</B//<<BBB//<0/<//<-<<CAC.C<.::;@@?-9-?9C??-@@B?-@-99--999-;-9-99-9---///-9----;--9--9-9---9////:A-;--9-A---9-9-/--@--/9///--@--9-99A-@-------//9-------////99/;/;9--9////;////;
############################## 
##############maybe need to make a foreach for all lanes
my $in= "/cs/cbio/oren/SequencingData/Hagar/Hagar3_S2_L001_R2_001.fastq";
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
my @anti=("GCC","GTT","GTC","GCT","TCC","TTT","TTC","TCT");
my @sense=("GCC","ACC","TCC");
foreach my $numrep (1..1){
	print "number of repeats is: ".$numrep."\n";
	my $out= "/cs/cbio/oren/SequencingData/Hagar/Hagar3.fa";
	my $out1= "/cs/cbio/oren/SequencingData/Hagar/Hagar3_anti".$numrep.".txt";
	my $out2= "/cs/cbio/oren/SequencingData/Hagar/Hagar3_sense".$numrep.".txt";
	my $out3= "/cs/cbio/oren/SequencingData/Hagar/Hagar3_oops".$numrep.".txt";
	open OUT,">$out";
	open OUT1,">$out1";
	open OUT2,">$out2";
	open OUT3,">$out3";
	my $oopscount=0;
	my $count=0;
	my @count;
	$count[0]=0;
	foreach my $n (1..($numrep-1)){
		$count[$n]=$count[$n-1]+3;
		#print $count[$n]."\n";
	}
	foreach my $key (keys(%hash)){
		print OUT ">".$key."\n".$hash{$key}."\n";
		#print $hash{$key}."\n";
		my $a=substr($hash{$key}, 83, ($numrep*3));
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
			$c[$x]=$b[$count[$x]].$b[$count[$x]+1].$b[$count[$x]+2];
			#print $c[$x]."\t";
			if ($c[$x] eq $anti[0]){$anticount++;}
			if (($c[$x] eq $anti[1])||($c[$x] eq $anti[2])||($c[$x] eq $anti[3])){$anticount++; $flagAnti=1;}
			if ($c[$x] eq $sense[0]){$sensecount++;}
			if (($c[$x] eq $sense[1])){$sensecount++; $flagSense=1;}	
		}
		#print $anticount."\t".$sensecount."\n";
		my $b=substr($hash{$key}, (83+($numrep*3)), 3);#GCG
		if(($b eq "ACA") || ($b eq "GCG") || ($b eq "GTG") || ($b eq "ACG") || ($b eq "GTA") || ($b eq "GCA")){
			#print $b."\n";
			if (($anticount==$numrep) && ($sensecount==$numrep)){print OUT3 ">".$key."\n".$hash{$key}."\n"; $count++;}
			if (($anticount==$numrep)&&($flagAnti==1)){$count++; print OUT1 ">".$key."\n".$hash{$key}."\n";}
			if (($sensecount==$numrep)&&($flagSense==1)){$count++; print OUT2 ">".$key."\n".$hash{$key}."\n";}
		}
		
	}
	close(OUT1); close(OUT2); close(OUT3); close(OUT);
	print "number: ".$count."\n";
	#print "oopscount: ".$oopscount."\n"."END\n";
}
