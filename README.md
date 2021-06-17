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


# Example

This example was designed to work on Euler, a bsub system

First setup a directory to work in.
```
mkdir Example
cd Example
mkdir inputDirectory
```
Now install/donwload the GBSmode scripts from github.

```
module load git
git clone https://github.com/stevenandrewyates/GBSmode
```
Next download some data. In this case 100,000 reads from 101 cassava samples will be downloaded. This is a slight hack as the downloaded data is directly summarised for GBSmode.
```
for x in $(seq 1717931 1718031); do echo "$HOME/sratoolkit.2.11.0-centos_linux64/bin/fastq-dump --split-files -X 100000 -Z SRR$x |grep -A1 \"^@\" |  grep -v "-" | grep -v "@" | cut -c-75 | sort | uniq -c > inputDirectory/$x.fastq";done
```
Now we need quite a bit of RAM. Because all the read numbers will be tabulated for the population. Being as this is a bsub system we need to clean up the ouptut. Please change the lsf job!
```
bsub -R "rusage[mem=4096]" -n 20 perl GBSmode/getUniqueTags.pl inputDirectory/ > Count.data
grep ^Tag lsf.o175615454  > Count.data
sed '1,/^Tag/d' lsf.o175615454  >> Count.data
```
Install the genome reference.
```
module load git
git clone https://github.com/stevenandrewyates/SAYReadMappingDNA
sh SAYReadMappingDNA/01_DownloadGenome.sh -f ftp://ftp.ensemblgenomes.org/pub/plants/release-50/fasta/manihot_esculenta/dna/Manihot_esculenta.Manihot_esculenta_v6.dna.toplevel.fa.gz
```
Read map using bowtie2.
```
bowtie2 -x GENOME/genome -U sequences.fastq -S file.sam
```
Filter the data, by selecting only linkage groups ("LG") and get the first six columns. The other columns just cause problems...
```
grep "LG" file.sam | grep -v "^@" | cut -f 1-6 | grep 'Hap' > filter.sam
```
Now remove unmapped data.
```
perl GBSmode/TrimCountsInput.pl filter.sam Count.data > filter.data
```
Use the GBSmode script to find polymorphisms.
```
R --vanilla --slave "--args filter.data filter.sam output fastq GENOME/genome.fasta" < GBSmode/GBSmode.R
```
For ease of use we will convert the GBSmode output to vcf format, so it is portable with other analysis programs 
```
R --vanilla --slave "--args output_genotype.txt output_genotype_coverage.txt GBSmode.vcf 1" <  GBSmode/ModeToVCF.R
```
To finish off we will prune the data and get loci with at least fifty genotype calls
```
git clone https://github.com/stevenandrewyates/filterVCF
R --vanilla --slave "--args GBSmode.vcf bowtieclean.vcf 5 50" < filterVCF/filterVCF.R
```
