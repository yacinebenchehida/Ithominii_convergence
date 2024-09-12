#! /usr/bin/perl -w
#Author: Prof. Kanchon Dasmahapatra
use strict;
use warnings;

my $group1 = $ARGV[0] or die "provide group1 as first option in this format indivN,indivN,indivN.\n";
my $group2 = $ARGV[1] or die "provide group2 as second option in this format indivN,indivN,indivN.\n";
my $min_ind1 = $ARGV[2] or die "provide minimum number of individual in group 1 as third option.\n";
my $min_ind2 = $ARGV[3] or die "provide minimum number of individual in group 2 as 4th option.\n";
my $input_file = $ARGV[4] or die "provide input vcf file as 5th option.\n";

my $output_calls = $input_file . '_minind' . $min_ind1 . '_' . $min_ind2 . '_calls'; 
my $output_fixed_diffs = $input_file . '_minind' . $min_ind1 . '_' . $min_ind2 . '_fixed_diff';
my $output_hom_het = $input_file . '_minind' . $min_ind1 . '_' . $min_ind2 . '_hom_het';

my $i; my $line; my @value; my @fields;
my @names; my $ntaxa = 0;

open INPUTFILE, "<$input_file" or die "could not open input file.\n";
open (OUTPUTCALLS, ">$output_calls");
open (OUTPUTFD, ">$output_fixed_diffs");
open (OUTPUTHOMHET, ">$output_hom_het");

my @group1_array = split(',',$group1);
my @group2_array = split(',',$group2);
my $group1_ntaxa = @group1_array;
my $group2_ntaxa = @group2_array; 

# Skip the header section of the vcf
while ($ntaxa == 0) {
	$line = <INPUTFILE>;
	chomp($line);
	@value = split(' ',$line);
	if ($value[0] eq '#CHROM') {
# Get number of samples in vcf file
		$ntaxa = @value - 9;
# Print the header line
		print OUTPUTCALLS "$value[0]\t$value[1]"; print OUTPUTFD "$value[0]\t$value[1]"; print OUTPUTHOMHET "$value[0]\t$value[1]";

		for ($i = 0; $i < $group1_ntaxa ; $i++) {
			$names[$i] = $value[8 + $group1_array[$i]];
			print OUTPUTCALLS "\t$names[$i]"; print OUTPUTFD "\t$names[$i]"; print OUTPUTHOMHET "\t$names[$i]";
		}
		print OUTPUTCALLS "\t"; print OUTPUTFD "\t"; print OUTPUTHOMHET "\t";
		for ($i = 0; $i < $group2_ntaxa ; $i++) {
			$names[$i] = $value[8 + $group2_array[$i]];	
			print OUTPUTCALLS "\t$names[$i]"; print OUTPUTFD "\t$names[$i]"; print OUTPUTHOMHET "\t$names[$i]";
		}
		print OUTPUTCALLS "\n"; print OUTPUTFD "\n"; print OUTPUTHOMHET "\n";
	}
}

while ($line = <INPUTFILE>) {
	chomp($line);    
	@value = split(' ',$line);
#	print OUTPUTCALLS "$value[0]\t$value[1]";
	my $group1_hom1 = 0; my $group1_hom2 = 0; my $group1_het = 0; my $group1_genos = 0; # setting the counters for group1
	my $group2_hom1 = 0; my $group2_hom2 = 0; my $group2_het = 0; my $group2_genos = 0; # setting the counters for group2
	my $null = 0;
	for ($i = 0; $i < $group1_ntaxa ; $i++) {                     # counting homozygotes and heterozygotes for group1
		@fields = split(':',$value[8 + $group1_array[$i]]);
		if ($fields[0] eq '0/0' || $fields[0] eq '0|0') {$group1_hom1++;}
		elsif ($fields[0] eq '1/1' || $fields[0] eq '1|1') {$group1_hom2++;}
		elsif ($fields[0] eq '0/1' || $fields[0] eq '0|1') {$group1_het++;}
		elsif ($fields[0] eq '1/0' || $fields[0] eq '1|0') {$group1_het++;}
		elsif ($fields[0] eq './.' || $fields[0] eq '.|.') {$null++;}
		else {print "unexpected genotype $fields[0] at $value[0] $value[1]\n";}
#		$fields[0] =~ s/\///g;
#		print OUTPUTCALLS "\t$fields[0]";
	}

	for ($i = 0; $i < $group2_ntaxa ; $i++) {                     # counting homozygotes and heterozygotes for group2
		@fields = split(':',$value[8 + $group2_array[$i]]);
                if ($fields[0] eq '0/0' || $fields[0] eq '0|0') {$group2_hom1++;}
                elsif ($fields[0] eq '1/1' || $fields[0] eq '1|1') {$group2_hom2++;}
                elsif ($fields[0] eq '0/1' || $fields[0] eq '0|1') {$group2_het++;}
                elsif ($fields[0] eq '1/0' || $fields[0] eq '1|0') {$group2_het++;}
                elsif ($fields[0] eq './.' || $fields[0] eq '.|.') {$null++;}
		else {print "unexpected genotype $fields[0] at $value[0] $value[1]\n";}
#		$fields[0] =~ s/\///g;
#		print OUTPUTCALLS "\t$fields[0]";
	}
#	print OUTPUTCALLS "\n";

	$group1_genos = $group1_hom1 + $group1_hom2 + $group1_het;
	$group2_genos = $group2_hom1 + $group2_hom2 + $group2_het;
	
	if (($group1_genos >= $min_ind1) && ($group2_genos >= $min_ind2)) {             # checking that both groups have at least $min_ind genotypes
		print OUTPUTCALLS "$value[0]\t$value[1]";
		for ($i = 0; $i < $group1_ntaxa ; $i++) {
			@fields = split(':',$value[8 + $group1_array[$i]]);
			$fields[0] =~ s/\///g;
                        $fields[0] =~ s/\|//g;
			print OUTPUTCALLS "\t$fields[0]";
		}
		for ($i = 0; $i < $group2_ntaxa ; $i++) {
                        @fields = split(':',$value[8 + $group2_array[$i]]);
                        $fields[0] =~ s/\///g;
                        $fields[0] =~ s/\|//g;
                        print OUTPUTCALLS "\t$fields[0]";
                }
		print OUTPUTCALLS "\n";

		# checking for positions that are fixed differences
		if ((($group1_hom1 == $group1_genos) && ($group2_hom2 == $group2_genos)) || (($group1_hom2 == $group1_genos) && ($group2_hom1 == $group2_genos))) {
			print OUTPUTFD "$value[0]\t$value[1]";
		        for ($i = 0; $i < $group1_ntaxa ; $i++) {                     # counting homozygotes and heterozygotes for group1
		                @fields = split(':',$value[8 + $group1_array[$i]]);
       		         	$fields[0] =~ s/\///g;
	                        $fields[0] =~ s/\|//g;
		                print OUTPUTFD "\t$fields[0]";
	        	}

		        for ($i = 0; $i < $group2_ntaxa ; $i++) {                     # counting homozygotes and heterozygotes for group2
        		        @fields = split(':',$value[8 + $group2_array[$i]]);
                		$fields[0] =~ s/\///g;
                                $fields[0] =~ s/\|//g;
                		print OUTPUTFD "\t$fields[0]";
        		}
	        	print OUTPUTFD "\n";
		} 

		# checking for hom het positions
		if (((($group1_hom1 + $group1_het) == $group1_genos) && ($group2_hom2 == $group2_genos)) || 
		    ((($group1_hom2 + $group1_het) == $group1_genos) && ($group2_hom1 == $group2_genos)) ||
		    (($group1_hom2 == $group1_genos) && (($group2_hom1 + $group2_het) == $group2_genos)) ||
		    (($group1_hom1 == $group1_genos) && (($group2_hom2 + $group2_het) == $group2_genos))) {
			print OUTPUTHOMHET "$value[0]\t$value[1]";
       	        	for ($i = 0; $i < $group1_ntaxa ; $i++) {                     # counting homozygotes and heterozygotes for group1
                        	@fields = split(':',$value[8 + $group1_array[$i]]);
	                        $fields[0] =~ s/\///g;
                                $fields[0] =~ s/\|//g;
	                        print OUTPUTHOMHET "\t$fields[0]";
	                }

	                for ($i = 0; $i < $group2_ntaxa ; $i++) {                     # counting homozygotes and heterozygotes for group2
	                        @fields = split(':',$value[8 + $group2_array[$i]]);
	                        $fields[0] =~ s/\///g;
                                $fields[0] =~ s/\|//g;
	                        print OUTPUTHOMHET "\t$fields[0]";
	                }
	                print OUTPUTHOMHET "\n";
		}
	}
}

close INPUTFILE; close OUTPUTCALLS; close OUTPUTFD; close OUTPUTHOMHET;
