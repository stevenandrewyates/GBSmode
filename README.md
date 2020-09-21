# GBSmode

## Description:

GBSmode is a pipeline for genotype calling using GBS data.

## Usage:

GBSmode starts with a folder with deplexed reads, called DEPLEX. Being as GBS barcode size varies, it trims the reads to a uniform length and the frequency of each unique sequence is counted. A series of ‘grep’ command filter the fastq file to retrieve the sequence data. The sequence data is then trimmed to a uniform length using the ‘cut’ command, in this case to a length of 75 bp. The output is then sent to the inputDirectory.

    mkdir inputDirectory
    for x in $(ls DEPLEX/); do grep -A1 "^@" DEPLEX/$x | grep -v "-" | grep -v "@" | cut -c-75 | sort | uniq -c > inputDirectory/$x;done

The sample read counts are summarised per population, generating a matrix of unique sequences (rows) with genotype counts (columns), using perl.

    perl getUniqueTags.pl inputDirectory/ > Count.data

the sequence data is then converted back to fastq format.

    perl makeFastq.pl Count.data > sequences.fastq

A genome reference is then prepared for read mapping

    bowtie2-build genome.fasta genome.ref

then the reads are mapped to the genome

    bowtie2 –x genome.ref –U sequences.fastq –S file.sam

the resulting ‘sam’ file is then filtered

    grep “Chr” file.sam | grep –v “^@” | cut –f 1-6 | grep ‘Hap’ > filter.sam

the count data is then filtered for reads that align to the genome, its’ arguments are the filtered ‘sam’ file and the count data.

    perl TrimCountsInput.pl filter.sam Count.data > filter.data

Genotyping is made using the “GBSmode.R’ program which can be run interactively in R or via a terminal. The arguments are the filtered data (filter.data), filtered sam file (filter.data), a prefix for the output files (output.prefix), a common file identifier (fastq), the genome reference used (genome.fasta).  For example:

    R –vanilla –slave “—args filter.data filter.sam output fastq genome.fasta” < GBSmode.R

The data can then be converted into fastq data using the following command. It takes the two output files from GBSmode.R (output_genotype.txt output_genotype_coverage.txt), an output file name (GBSmode.vcf) and the number of cores to use (1 in this example).

    R --vanilla --slave "--args output_genotype.txt output_genotype_coverage.txt GBSmode.vcf 1" <  ModeToVCF.R
