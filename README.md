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

# Easy example

This example is designed for a quick check and uses only the dependencies needed for GBSmode

1) Download the data
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/242/695/GCF_000242695.1_LepOcu1/GCF_000242695.1_LepOcu1_genomic.fna.gz
gunzip GCF_000242695.1_LepOcu1_genomic.fna.gz
```

2) Prepare the genome for bowtie2-build
```
mkdir GENOME
```

3) Extract the first two linkage groups, less time later
```
head -n 1600000 GCF_000242695.1_LepOcu1_genomic.fna | sed "s/>.*group />/g" > GENOME/genome.fasta
bowtie2-build GENOME/genome.fasta GENOME/genome
```

4) Download the raw data
```
mkdir DEPLEX
cd DEPLEX
for x in $(seq 1 20); do echo wget https://zenodo.org/record/1219888/files/progeny_$x;done | sh
cd ..
```

5) Trim the data to 70 base pairs and then sort and count the frequency of each sequence
prepare and sort the data
```
mkdir inputDirectory
for x in $(ls DEPLEX/); do grep -A1 "^>" DEPLEX/$x | grep -v "-" | grep -v ">" | cut -c-70 | sort | uniq -c > inputDirectory/$x;done
```

6) Download GBSmode 
```
git clone https://github.com/stevenandrewyates/GBSmode
```

7) Prepare data for GBSmode
```
perl GBSmode/getUniqueTags.pl inputDirectory/ > Count.data
perl GBSmode/makeFastq.pl Count.data > sequences.fastq
bowtie2 -x GENOME/genome -U sequences.fastq -S file.sam
grep "LG" file.sam | grep -v "^@" | cut -f 1-6 | grep 'Hap' > filter.sam
perl GBSmode/TrimCountsInput.pl filter.sam Count.data > filter.data
```
8) Lower the minimum number of genotypes to 10
```
sed -i 's/NG <- 19/NG <- 9/' GBSmode/GBSmode.R
```
10) Run GBSmode
```
R --vanilla --slave "--args filter.data filter.sam output progeny GENOME/genome.fasta" < GBSmode/GBSmode.R
```
11) Convert the GBSmode data into vcf format 
```
R --vanilla --slave "--args output_genotype.txt output_genotype_coverage.txt GBSmode.vcf 1" <  GBSmode/ModeToVCF.R
```
12) Voila, the output VCF is: "GBSmode.vcf"


# Example

This example was designed to work on Euler, a bsub system

First setup a directory to work in.
```
mkdir Example
cd Example
mkdir inputDirectory
```
Now install/download the GBSmode scripts from github.

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
