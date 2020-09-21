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

# Perl script for filtering the count data based on
# reads which align

# its’ arguments are the filtered ‘sam’ file and 
# the count data

# output is written to standard output,
# so needs piping to file

##################################################
######              Usage                   ######
##################################################

# perl TrimCountsInput.pl filter.sam Count.data > filter.data

##################################################
######              Script                  ######
##################################################



my @files;	#also sample IDs
my %haps;
my $counter = 0;

### read in the SAM file 

open(INPUT, @ARGV[0]) or die "could not open file";
while (<INPUT>) 
{
chomp;
my $line =$_;
my @data = split("\t",$line);
#print $data[0];
$data[0] =~ s/@//;
$haps{$data[0]} = 0;
}
close(INPUT);

# now read the counts file
open(COUNT, @ARGV[1]) or die "could not open file";
while (<COUNT>) 
{

chomp;
my $line =$_;
if ($counter == 0) {print $line, "\n"};
$counter++;

my @data = split("\t",$line);
$data[0] =~ s/@//;
#print $data[0],"\n";
if (exists($haps{$data[0]})) { print $line,"\n";}
}
close(COUNT);

