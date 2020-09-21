
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

# Script for genotyping GBS data

# it needs five 'args' listed consecutively below
# 1) the count data 
# 2) a filtered sam file
# 3) output file name prefix
# 4) a delimeter for the sample names (like fastq)
# 5) the reference genome, in fasta format

# you can also specify biological factors in the 
# population here.

# ploidy state of the genotypes
ploidy <- 2

# the maximum number of haplotypes within the population,
# Four assumes a bi-parental population with heterozygous parents
maxhaps <- 4

# number of genotypes, to be considered before counting a mutation
NG <- 19 

# minimum allele frequenc 
# min read depth for an allele at a loci
minAF <- 100

# minimum allele count, 
# the number of unique reads before an individual is considered homozygous
minAC <- 8

# the script will write data to a log file throughout
# it will produce two output files,
# one with the genotyping
# another with the corresponding count data.

# it will also produce some useful graphs
# which can be used to see the coverage per loci
# this might be useful to see how well the library 
# worked :)

##################################################
######              Usage                   ######
##################################################

# the script can be run interactively, in R or via
# command line:

# R -vanilla -slave "-args filter.data filter.sam output fastq genome.fasta" < GBSmode.R

##################################################
######             Requirements             ######
##################################################


library("reshape2")
library("ggplot2")
library("seqinr")
library("GenomicAlignments")


##################################################
######              Script                  ######
##################################################

# input requirements
args <- commandArgs(trailingOnly = T)


COUNTFILE <- args[1]
SAMFILE <- args[2]
OUTPUTFILE <- args[3]
FILEDELIMITER <- args[4]
GENOME <- args[5]

#################################
# give some information
logfile <- paste(OUTPUTFILE,".log.txt",sep="")

write.table(file=logfile,print(paste("log file is = ",logfile,sep="")),col.names=F,row.names=F,quote=F)
print(paste("log file is = ",logfile,sep=""))

print(paste("Count file =",COUNTFILE))
write.table(file=logfile,print(paste("Count file =",COUNTFILE)),col.names=F,row.names=F,quote=F,append=T)

print(paste("SAM file =",SAMFILE))
write.table(file=logfile,print(paste("SAM file =",SAMFILE)),col.names=F,row.names=F,quote=F,append=T)

print(paste("Output file =",OUTPUTFILE))
write.table(file=logfile,print(paste("Output file =",OUTPUTFILE)),col.names=F,row.names=F,quote=F,append=T)


print(paste("File delimter =",FILEDELIMITER))
write.table(file=logfile,print(paste("File delimter =",FILEDELIMITER)),col.names=F,row.names=F,quote=F,append=T)


#########################################
# write the biological factors to the log

write.table(file=logfile,print(paste("ploidy = ",ploidy)),col.names=F,row.names=F,quote=F,append=T)

write.table(file=logfile,print(paste("Maximum number of haplotypes = ",maxhaps)),col.names=F,row.names=F,quote=F,append=T)

write.table(file=logfile,print(paste("minimum number of genotypes = ",NG+1)),col.names=F,row.names=F,quote=F,append=T)

write.table(file=logfile,print(paste("minimum allele frequency = ",minAF)),col.names=F,row.names=F,quote=F,append=T)

write.table(file=logfile,print(paste("minimum allele count = ", minAC)),col.names=F,row.names=F,quote=F,append=T)

#################################
# read in the data


# read in the genome
write.table(file=logfile,print(paste("reading in",GENOME)),col.names=F,row.names=F,quote=F,append=T)
genome <- read.fasta(file=GENOME)
write.table(file=logfile,print(paste("read in",GENOME)),col.names=F,row.names=F,quote=F,append=T)

# read the data
c <- read.table(COUNTFILE,sep="\t",header=T)
print(paste("Length of Count file =",length(c[,1])))
write.table(file=logfile,print(paste("Length of Count file =",length(c[,1]))),col.names=F,row.names=F,quote=F,append=T)

# s <- read.table("raw1.sort.sam",sep="\t")
s <- read.table(SAMFILE,sep="\t")
print(paste("Length of SAM file =",length(s[,1])))
write.table(file=logfile,print(paste("Length of SAM file =",length(s[,1]))),col.names=F,row.names=F,quote=F,append=T)

# rename the columns
colnames(s)[1:6] <- c("ID","Flag","Ref","Pos","MAPQ","CIGAR")
si <- s[,c(1:6)]

# format the rownames
c$ID <- gsub('@','',row.names(c))
si$ID <- as.character(si$ID)

# get the correct starting position
sf <- cbind(si, as.data.frame(unlist(extractAlignmentRangesOnReference(s$CIGAR, pos=s$Pos))))
si$Pos[sf$Flag == 16] <- sf$end[sf$Flag == 16] 

# merge the datasets
da <- merge(si,c,by="ID")

################################
# just check what happened

print(paste("Length after merging SAM and Count files =",length(da[,1])))
write.table(file=logfile,print(paste("Length after merging SAM and Count files =",length(da[,1]))),col.names=F,row.names=F,quote=F,append=T)

# find duplicated read mappings

dups <- da$ID[duplicated(da$ID)]
print(paste("Number of Duplicates found =",length(dups)))
write.table(file=logfile,print(paste("Number of Duplicates found =",length(dups))),col.names=F,row.names=F,quote=F,append=T)

# remove the duplicates
dc <- da[!(da$ID %in% dups),]
print("Removing the duplicates")
write.table(file=logfile,print("Removing the duplicates"),col.names=F,row.names=F,quote=F,append=T)

# find the unique positions
dc$Upos <- paste(dc$Flag,dc$Ref,dc$Pos,sep="_")

# make a vector of this
UPOS <- unique(dc$Upos)
print(paste("Number of unique tag origins =",length(UPOS)))
write.table(file=logfile,print(paste("Number of unique tag origins =",length(UPOS))),col.names=F,row.names=F,quote=F,append=T)

####################################
# get a list of useful positions

over100 <- dc$Upos[dc$Total > minAF]
over100 <- unique(over100)
print(paste("Number of alleles with over", minAF, "reads =",length(over100)))
write.table(file=logfile,print(paste("Number of alleles with over", minAF, "reads =",length(over100))),col.names=F,row.names=F,quote=F,append=T)

# set up something to store the data
OUT <- NULL
OUT2 <- NULL

###################################################################
# start to call the genotypes
print("Calling the genotypes this may take some time")
write.table(file=logfile,print("Calling the genotypes this may take some time"),col.names=F,row.names=F,quote=F,append=T)

# used for counting to output the progress
dummycount <- 0;

##############################################
##############################################
# start of loop
##############################################
##############################################

for (p in 1:length(over100))
{
# start by incrementing the locis processed
dummycount <- dummycount +1
# if the number is divisible by 1000 then output this information
if ((dummycount %% 1000) == 0) 
{
	print(paste("Processed",dummycount,"of",length(over100),"loci, found",dim(OUT)[1],"informative tags"))
	write.table(file=logfile,print(paste("Processed",dummycount,"of",length(over100),"loci, found",dim(OUT)[1],"informative tags")),col.names=F,row.names=F,quote=F,append=T)
}

# select a loci which is useful
l <- dc[dc$Upos == over100[p],]
l$Tag <- as.character(l$Tag)

################################################
# get the informative haplotypes, with 5% minor allele frequency count 
Itags <- unique(as.character(l$Tag[l$Total > ((sum(l$Total)/100)*5)]))

# split these into individual bases
if (length(Itags) < 2) 
{next}
ICtags <- strsplit(Itags,"")

############################################
# if it isn't polymorphic then skip, i.e. only one haplotype
if (length(Itags) == 1) 
{next}

############################################
# if there are more than ploidy (2) then check it isn't a from a multimap
if (length(Itags) > ploidy) 
	{
	# make a quick scan to see if a genotype has more than the ploidy level of haplotypes
	zu <- l[l$Tag %in% Itags,grep(FILEDELIMITER,colnames(l))]
	# make a vector to store the data
	QuickHap <- NULL
	# get the maximum number of observed haplotypes
	for (u in 1:length(zu)) {QuickHap <- c(QuickHap,sum(zu[,u] > 0))}
	# check if any are greater than the ploidy level, if it is then skip
	if (any(QuickHap > ploidy)) {next}
	}

############################################
# check that the number of alleles does not exceed the number expected haplotypes
if (length(Itags) > maxhaps) {next}

###################################################
# get the genomic sequence
GP <- l$Pos[1]
GN <- names(genome)[l$Ref[1] == names(genome)]
GI <- max(nchar(as.character(l$Tag)))-1
GM <- paste(max(nchar(as.character(l$Tag))),"M",sep="")
GS <- toupper( paste(genome[[GN]][GP:(GP+GI)],collapse=""))

######################################################
# switch the strting position to the origin of the read if it is reversed (Flag 16)
if(l$Flag[1] == 16)
{
if(GP-GI < 1) {next} 
GS <- toupper(paste(comp( paste(genome[[GN]][GP:(GP-GI)])),collapse=""))
}

# prepare for alignement
qseq <- DNAStringSet(c(as.character(l$Tag),GS))
cigar <- c(as.character(l$CIGAR),GM)

# if the read is reversed switch the cigar string
if(l$Flag[1] == 16)
{
for(CCF in 1:length(cigar))
{
SE <- strsplit(cigar,"[A-Z]")[[CCF]]
SN <- strsplit(cigar,"[0-9]")[[CCF]]
SN <- SN[nchar(SN) >0]
cigar[CCF] <- gsub(" ","",paste(rev(SE),rev(SN),collapse=""))
}
}

########################################################
# make the pairwsie alignment INSERTIONS will be omitted

clipped_qseq <- sequenceLayer(qseq, cigar,from="query", to="query-after-soft-clipping")
WE <- sequenceLayer(clipped_qseq, cigar, from="query", to="reference")

# get the alignment as character
l$ALN <- as.character(WE[1:(dim(l)[1])])
# trim them all to the same (smallest) size
TL <-  min(nchar(l$ALN[1:(length(l$ALN))]))
l$ALN <- substr(l$ALN,0,TL)
# and trim the cigar
cigar <- cigarNarrow(cigar[1:(length(l$ALN))], end=TL)[1:(dim(l)[1])]


#genomic sequence
GA <- as.character(WE[length(WE)])

# now just get the genotype data
gd <- l[,grep(FILEDELIMITER,colnames(l))]

# variable to store potential alleles
informative <- NULL

######################################################################	
# get the informative haplotypes, with 5% minor allele frequency count 

l$CI <- cigar
IOUT <- l[l$Total > ((sum(l$Total)/100)*5),c("Tag","ALN","CIGAR","CI")]
# concatenate the aligned sequence and the condensed cigar
IOUT$COM <- paste(IOUT$ALN,IOUT$CI)
# now remove duplications
IOUT <- IOUT[!duplicated(IOUT$COM) ,]
if(dim(IOUT)[1] < 2) {next}
Itags <- IOUT$ALN

# split these into individual bases
ICtags <- strsplit(Itags,"")
Icigar <- IOUT$CI

Isite <- NULL
u <- 1:nchar(l$ALN[1])
for (t in 2:length(ICtags)) {Isite <- c(Isite,u[(ICtags[[1]] == ICtags[[t]]) == FALSE])}
Isite <- unique(Isite)
Isite <- Isite[order(Isite)]
if(length(Isite) < 1 & length(unique(Icigar)) < 1) {next}
UHAP <- NULL
IOUT$P <- NULL
for (c in 1:length(ICtags)) 
	{
	UHAP <- c(UHAP, paste(paste(paste(Isite,ICtags[[c]][Isite],sep=""),collapse=""),Icigar[c],sep=""))
	IOUT$P[c] <- paste(paste(Isite,ICtags[[c]][Isite],sep=""),collapse="")
	}

############################################

############################################
# if it isn't polymorphic then skip, i.e. only one haplotype
if (length(UHAP) == 1) {next}

# or it has more haplotypes than possible in the population
if (length(UHAP) > maxhaps) {next}

#############################################
# if it isn't call the the alleles
if (length(Itags) > 1)
	{
	# make the original tag characters
	l$ALN<- as.character(l$ALN)

	# make a simple count
	u <- 1:nchar(l$ALN[1])
	
	# now make the 
	l$COM <- NULL
	ICtags2 <- strsplit(l$ALN,"")
	
	# extract the polymorphic positions
	for (t in 1:length(l$ALN)) {l$COM[t] <- paste(paste(Isite,strsplit(l$ALN[t],"")[[1]][Isite],sep=""),collapse="")}
	
	# combine with the reduced cigar
	l$COM2 <- paste(l$COM,l$CI,sep="")

	# check that it is polymorphic
	if(length(Isite) < 1 & length(unique(Icigar)) < 1) {next}

	# now remove the tags which did not match the expected haplotypes
	l <- l[l$COM2 %in% UHAP,]

	# now get the haplotype frequencies
	COUNTS <- colSums(l[l$COM2 == UHAP[1],grep(FILEDELIMITER,colnames(l))])
	
	# loop through the remaining haplotypes
	for (d in 2:length(UHAP)) {COUNTS <- rbind(COUNTS, colSums(l[l$COM2 == UHAP[d],grep(FILEDELIMITER,colnames(l))]))}

	#############################################
	# record how many loci were scorable
	
	# instantiate a variable for the output
	HA <- NULL
	
	# for recording the counts
	CCT <- paste(COUNTS)
	CRS <- NULL
	for (w in 1:(dim(COUNTS)[2])) {CRS <- c(CRS,paste(COUNTS[,w],collapse="_"))}
	SCF <- IOUT$P
	
	# prepare the locus data for output (with the other loci)
	# transpose the data (wide) and add the locus information
	KI <- 1:length(Icigar)
	KI[grep("I",Icigar)]

	# add indels
	for(j in 	KI[grep("I",Icigar)])
	{
	VARPOS <- data.frame(CS=cigarRangesAlongQuerySpace(Icigar, with.ops=TRUE)[[j]][,1]@start)
	VARPOS$CW <- cigarRangesAlongQuerySpace(Icigar, with.ops=TRUE)[[j]][,1]@width
	VARPOS$CL <- (VARPOS$CS + VARPOS$CW) -1
	VARPOS$TY <- cigarRangesAlongQuerySpace(Icigar, with.ops=TRUE)[[j]][,1]@NAMES
	VARPOS$IN <- NA
	VARPOS <- VARPOS[VARPOS$TY == "I",]
	for (g in 1:length(VARPOS[,1])) {VARPOS$IN[g] <- substr(IOUT$Tag[j],VARPOS$CS[g],VARPOS$CL[g])}
	VARPOS$GP <- VARPOS$CL-cumsum(VARPOS$CW)
	VARPOS$INDEL <- paste(VARPOS$GP,VARPOS$TY,VARPOS$IN,sep="")
	SCF[j] <-  paste(SCF[j],paste(VARPOS$INDEL,collapse=""),sep="")
	print("added indel")
	print(p)
	}
	SCF2 <- SCF
	SCF <- paste(SCF,collapse="_")

	ko2 <- as.matrix(t(c(over100[p],paste(Isite[order(Isite)],collapse="_"),SCF,CRS,as.character(l$Tag[l$Total == max(l$Total)])[1],GA)))
	
	# add column names for formatting
	ko2 <- as.data.frame(ko2)
	colnames(ko2) <- c("LOCI","Isite","Alleles",colnames(l)[grep(FILEDELIMITER,colnames(l))],"REFseq","GENseq")
	
	# change to a dataframe, for formatting
	ko2 <- as.data.frame(ko2)
	
	# add this locus to the other loci output
	OUT2 <- rbind(OUT2,ko2) 
	
	# voila

	###############################################
	
	# now call the alleles
	for (x in 1:(dim(COUNTS)[2]))
		{
		# set the non informative
		ha <- NA 

		HAPNUM <- 1:length(Itags)
		COM <- combn(HAPNUM,2)

		#############################################
		# first call the homozygotes
		for (k in 1:length(Itags)) 
			{
			if (COUNTS[k,x] >= minAC & sum(COUNTS[-k,x]) == 0) 
				{ha <-  paste(SCF2[k],SCF2[k],sep="_")}
			}

		##############################################
		# now call the heterozgotes
		for (k in 1:(dim(COM)[2])) 
			{
			if (COUNTS[(COM[1,k]),x] > 0 & COUNTS[(COM[2,k]),x] > 0) 
				{
				ha <-  paste(SCF2[[(COM[1,k])]],SCF2[[(COM[2,k])]],sep="_")
				}
			}
		
		# join the data
		HA <- c(HA,ha)
		}

	############################################
	# prepare the locus data for output (with the other loci)
	# transpose the data (wide) and add the locus information
	ko <- as.matrix(t(c(over100[p],paste(Isite[order(Isite)],collapse="_"),HA,as.character(l$Tag[l$Total == max(l$Total)])[1],GA)))

	# add column names for formatting
	colnames(ko) <- c("LOCI","Isite",colnames(l)[grep(FILEDELIMITER,colnames(l))],"REFseq","GENseq")
	
	# change to a dataframe, for formatting
	ko <- as.data.frame(ko)
	print(ko)
	
	# add this locus to the other loci output
	OUT <- rbind(OUT,ko) 
	
	# voila
	}
}

#####################################################
####   get some data on the gentoypes and loci   ####
#####################################################

# calculate how many genotypes are called per loci
FOUND <- NULL

# loop across the rows to do this
for (d in 1:length(OUT[,1]))
{

  # get a count of how many are not NA 
  FOUND <- c(FOUND,length(OUT[d,!is.na(OUT[d, grep(FILEDELIMITER,colnames(OUT))])]))
}

# convert the output to a dataframe, for plotting
FC <- data.frame(x=FOUND)

# put the information in the log file
write.table(file=logfile,print(paste("Writing Genotypes per loci data to ",OUTPUTFILE,"_GenotypesPerLoci.pdf",sep="")),col.names=F,row.names=F,quote=F,append=T)

# make a graph of the data
pdf(paste(OUTPUTFILE,"_GenotypesPerLoci.pdf",sep=""),width=6,height=6)
ggplot(FC, aes(x=x)) + geom_histogram(bins=20,colour="black", fill="grey")+theme_bw() +xlab("Number of genotypes per loci")
dev.off()

#################################################
# calculate how many loci are called per genotype
GF <- NULL

# loop across the columns
# get a count of how many are not NA 
for (d in 3:length(OUT[1,]))
{

  # get a count of how many are not NA 
  GF <- c(GF,length(OUT[!is.na(OUT[,d]),d]))
}

# convert the output to a dataframe, for plotting
GFC <- data.frame(x=GF)

# put the information in the log file
write.table(file=logfile,print(paste("Writing loci per genotype data to = ",OUTPUTFILE,"_LociPerGenotype.pdf",sep="")),col.names=F,row.names=F,quote=F,append=T)

# make a graph of the data
pdf(paste(OUTPUTFILE,"_LociPerGenotype.pdf",sep=""),width=6,height=6)
ggplot(GFC, aes(x=x)) + geom_histogram(bins=20,colour="black", fill="grey")+theme_bw() +xlab("Number of loci per genotype")
dev.off()

################################################
######    prepare the data for output  #########
################################################

# transpose the data
dr <- t(OUT)

# find the duplicates
dupx <- duplicated(dr[1,])

# remove duplicated
dx <- dr

# make column names
colnames(dx) <- paste(dx[1,],dx[2,],sep="_")
dx <- dx[c(-1,-2),]
dr <- dx
dr <- as.data.frame(dr)

SEQDATA <- as.data.frame(t(dr[(dim(dr)[1]),]))
SEQDATA$LOCI <- rownames(SEQDATA)
write.table(SEQDATA,file="seqdata.txt",sep="\t",quote=F,row.names=F)


dr$Original <- rownames(dr)

write.table(dr,file=paste(OUTPUTFILE,"_genotype.txt",sep=""),sep="\t",quote=F,row.names=F)

##########################

write.table(OUT2,file=paste(OUTPUTFILE,"_genotype_coverage.txt",sep=""),sep="\t",quote=F,row.names=F)


#####################################
# calculate the coverage by loci
wert <- data.frame(upos=UPOS,val=0)
for (t in 1:length(wert[,1])) {wert[t,2] <- sum(dc$Total[dc$Upos == wert[t,1]])}
plot(cumsum(wert$val[order(wert$val,decreasing=T)]))

COVERAGE <- data.frame(Nreads=cumsum(wert$val[order(wert$val,decreasing=T)]),Nloci=1:length(UPOS))

Q1 <- max(COVERAGE$Nloci[COVERAGE$Nreads < max(COVERAGE$Nreads)/4])
Q2 <- max(COVERAGE$Nloci[COVERAGE$Nreads < max(COVERAGE$Nreads)/2])
Q3 <- max(COVERAGE$Nloci[COVERAGE$Nreads < max(COVERAGE$Nreads)*0.75])


CHOP <- seq(1,length(COVERAGE[,1]),by=100)
COVERAGE <- COVERAGE[CHOP,]

pdf(paste(OUTPUTFILE,"_coverage.pdf",sep=""),width=6,height=6)
ggplot()+theme_bw()+ 
geom_rect(aes(xmin=0,ymin=0,xmax=Q3,ymax=(max(COVERAGE$Nreads)*0.75)),alpha=0.2,fill="red")+
geom_rect(aes(xmin=0,ymin=0,xmax=Q2,ymax=(max(COVERAGE$Nreads)*0.5)),alpha=0.2,fill="red")+
geom_rect(aes(xmin=0,ymin=0,xmax=Q1,ymax=(max(COVERAGE$Nreads)*0.25)),alpha=0.2,fill="red")+
geom_line(data=COVERAGE,aes(x=Nloci,y=Nreads))+ 
xlab(paste("Number of Loci\n(total ",max(COVERAGE$Nloci),")",sep=""))+
ylab(paste("Number of Reads\n(total ",max(COVERAGE$Nreads),")",sep=""))+
annotate("text",x=Q3*2,y=(max(COVERAGE$Nreads)*0.75),label=paste("Q3 =",Q3))+
annotate("text",x=Q3*2,y=(max(COVERAGE$Nreads)*0.5),label=paste("Q2 =",Q2))+
annotate("text",x=Q3*2,y=(max(COVERAGE$Nreads)*0.25),label=paste("Q1 =",Q1))+
ggtitle(paste(OUTPUTFILE,"_coverage.pdf"))
dev.off()


