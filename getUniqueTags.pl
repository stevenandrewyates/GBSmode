#!/usr/bin/perl
use strict;

##################################################
######            Details                   ######
##################################################

# author:   Steven Yates
# contact:  steven.yates@usys.ethz.ch
# Year:     2020
# Citation: TBA

##################################################
######            Description               ######
##################################################

# Perl script for converting a individual count 
# files to a population count file

##################################################
######              Usage                   ######
##################################################

# perl getUniqueTags.pl inputdirectory > output

##################################################
######              Script                  ######
##################################################

# set up variables for data collection
my %barcodes;
my @geno;

# read the directory from arguments

my $dirname=@ARGV[0];

# open the directory
opendir(DIR, $dirname);

# make an array of files in the directory
my @files = readdir(DIR);

# close the directory
closedir DIR;

# instantiate a counter variable at 0
my $Count = 0;

# for loop to iterate through the files
foreach my $key (@files)
 {

# add the file to the geno array
	push @geno, $key;

# create a variable with directory and file name 
	my $FI = $dirname."/".$key;

# open the file
	open(INPUT, $FI) or die "could not open file";

# read thorough the file
	while (<INPUT>) 
		{
# chomp the line
		chomp;
# assign the line to the line variable
		my $line =$_;
# split the data in the line by column (spaces)
		my @data = split(" ",$line);
# assign the barcode with the genotype and the count data
		$barcodes{$data[1]}{$key} = $data[0];
# update the total count for a barcode
		if (exists $barcodes{$data[1]}{Total}) {$barcodes{$data[1]}{Total} += $data[0]} else {$barcodes{$data[1]}{Total} = $data[0]};
# close the loop
		}
# close the file
	close($FI)
 } 


# print the header
print "Tag\t";

	foreach my $key (@files)
		{
#		if ($key =~ /GBS/) 
#			{
			print $key,"\t";
			
#			}
		}
	print "Total\n";  	



foreach my $B (sort keys %barcodes) 
	{
	if ($barcodes{$B}{Total} < 3) {next};#
	print "\@Hap_${Count}x\t",$B,"\t";
	$Count++;
	foreach my $key (@files)
		{
#		if ($key =~ /GBS/) 
#			{
			if (exists $barcodes{$B}{$key}) {print $barcodes{$B}{$key}, "\t"}  else {print "0\t"}
			
#			}
		}
	print $barcodes{$B}{Total};
	print "\n";  	
	}
