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

# Perl script for converting a count file to a fastq file

# it takes the first column as fastq header
# and the second as the sequence
# qualtiy scores are automatically assigned value of 40

# output is written to standard output,
# so needs piping to file

##################################################
######              Usage                   ######
##################################################

# perl makeFastq.pl input > output

##################################################
######              Script                  ######
##################################################

my @files;	#also sample IDs
my $counter = 0;

### read in the files for counting

open(INPUT, @ARGV[0]) or die "could not open file";
while (<INPUT>) 
{
$counter ++;
chomp;
my $line =$_;
if ($line =~ "^Tag") {next};
my @data = split("\t",$line);
#print $data[0],"\t",$data[1],"\t",$data[2],"\n";
#print "\@hap_${counter}_number=",$data[0],"\n",$data[1],"\n","+\n";
print $data[0],"\n",$data[1], "\n+\n";
#print "\@hap_${counter}x\n",$line,"\n","+\n";
my $L = length $data[1];
print "f"x $L;
print"\n";
}

