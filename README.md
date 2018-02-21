# Bisulfite_Footprinting
"Find_Repeats_FRM1_FASTQ.pl": will extract reads that contain 3bp repeats of the CGG and all possible conversion combinations. 
The repeats number is dictated by the first foreach loop.
The script will define "sense" vs "antisense" slipouts and all undefined/unmodified reads as "oops".
"Find_Repeats_C9_FASTQ.pl" is similar to "Find_Repeats_FRM1_FASTQ.pl" but for GGGGCC 6bp repeats conversion combinations.
"Find_RPL13A_FASTQ.pl" designed to work on control locus (non-repetative R-loop forming locus) to define "sense" vs "antisense" slipouts and undefined/unmodified reads as "oops".
Nucleotide_to_Number_based_coding.pl will produce a number based matrix according to conversion patterns relativle to the uncoverted reference sequence: 
C -> C = "1"
C -> T = "2"
C -> A/G = "9"
G -> G = "4"
G -> A = "5"
G-> T/C = "7"
In positions of which the refernce sequence is not G or C it will produce "0".
