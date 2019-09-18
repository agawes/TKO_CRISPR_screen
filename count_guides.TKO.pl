#!/usr/bin/perl
use strict;

my $fastq = $ARGV[0];

### read in the guide sequences
my $guides='/well/mccarthy/production/genetic-screens/data/EndoC_CRISPR_TKO/reference/tkov3_guide_sequence.txt';

my %guides;
open(IN,'<',$guides) or die;
my $head=<IN>;
while (my $line = <IN>){
	chomp $line;
	chop($line) if ($line =~ m/\r$/);

	my @line=split("\t",$line);
	$guides{$line[1]} = $line[0]."_".$line[2]."_".$line[3];
	#print "$line[1]\t$guides{$line[1]}\n";
}
close IN;

### read in the seq reads - use just R1
open(IN, "gunzip -c $fastq |") or die "gunzip $fastq: $!";
while (my $line = <IN>){ ## read header
	my $seq_line = <IN>;	
	### extract guide sequence from the read by matching for sequences of the U6 promoter and gRNA backbone
	### like in Extraction Pattern image at: https://cran.r-project.org/web/packages/caRpools/vignettes/CaRpools.html
	###  CACCG guide GTTTTAGAGC
	if ($seq_line =~ m/CACCG(.{20})GTTTTAGAGC/){
		print $1,"\t", $guides{$1},"\t";
		my @g=split("_",$guides{$1});
		print $g[0],"\n";
		
	} #else {print STDERR "Guide not found: $seq_line";}
	
	my $third_line=<IN>;
	my $q_line=<IN>;	
}
close IN;







## could extend it to also check if R2 has the reverse complement of this sequence, 
## but this is likely an overkill

### match the extracted guide sequence to the references