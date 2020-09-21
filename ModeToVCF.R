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

# Script for converting GBSmode output to vcf format.
# It takes the two output files from GBSmode.R 
# the genotype file (output_genotype.txt)
# the count data (output_genotype_coverage.txt), 
# an output file name (GBSmode.vcf) 
# and the number of cores to use (1 in this example).


##################################################
######              Usage                   ######
##################################################

# the script can be run interactively, in R or via
# command line:

# R --vanilla --slave "--args output_genotype.txt output_genotype_coverage.txt GBSmode.vcf 1" <  ModeToVCF.R

##################################################
######             Requirements             ######
##################################################

# if using multiple cores then:
# library(foreach)
# is required

##################################################
######              Script                  ######
##################################################

args <- commandArgs(trailingOnly = T)

d <- read.table(args[1],header=T,sep="\t")
COUNTDATA <- read.table(args[2],header=T,sep="\t")

CORE <- as.numeric(as.character(args[4]))

AC <- function(u) {as.numeric(as.character(u))}


makeVCF <- function(x)
{
	# CLEAR PREVIOUS OUTPUT
	VCF <- NULL
	SNPS <- NULL
	SLOC <- NULL
	# get the count data
	CR <- COUNTDATA[x,]
	# unlist the counts
	CRindividual <- unlist(strsplit(as.character(t(CR[4:(length(CR)-2)])),"_"))
	# get the count calls
	CD <- strsplit(as.character(CR[3][[1]]),"_")[[1]]
	# add missing insertion
	if(any(grep("_$",as.character(CR[3][[1]])))) { CD <- c(CD,"")}
	if(any(grep("^_",as.character(CR[3][[1]])))) { CD <- c(CD,"")}
	CD <- unique(CD)
	CF <- 1:2
		#length of data
		CL <- length(CR)
		REFseq <- as.character(t(CR[CL-1]))
		GENseq <- as.character(t(CR[CL]))
		# GET THE INFORMATION ON THE LOCUS
		LOC <- colnames(d)[x]
		# SPLIT THE DATA
		LOCS <- strsplit(LOC,"_")[[1]]
		# FIND THE DIRECTION OF THE READS (ORIENTATION)
		DIR <- LOCS[1]
		# GET THE SCAFFOLD/CHROMOSOME
		CHR <- LOCS[2]
		# GET THE OFFSET OF THE READ POSITION
		OFS <- LOCS[3]
		# GET THE LENGTH OF THE READ
		SIZ <- nchar(as.character(d[(dim(d)[1]),x]))
		SNPpos <- as.numeric(as.character(LOCS[4:length(LOCS)]))
		########################################
		# NOW PROCESS THE DATA DEPENDING UPON THE READ DIRECTION
		#######################################
		# negative offset
		if (DIR == "X16")
		{
		# FIND THE POSITIONS OF THE SNPS
		SLOC <-  AC(OFS)-SNPpos+1
		}
		# positive offset
		if (DIR == "X0")
		{
		# FIND THE POSITIONS OF THE SNPS
		SLOC <-  AC(OFS) -1 + SNPpos
		}


		##################################
		# get the reference SNPs
		REF <- NULL
		for (t in 1:length(SLOC)){ REF <- c(REF,substr(GENseq,SNPpos[t],SNPpos[t]))}

			##################################
			# FIND THE SNPS IN THE POPULATION
			# EXTRACT ALL UNIQUE HAPLOTYPES
			GE <- as.character(na.omit(d[1:(length(d[,x])-2),x]))
			# SPLIT THE OUTPUT
			GS <- strsplit(GE,"_")
			# NOW MAKE A UNIQUE LIST OF HAPLOTYPES
			GS <- unique(unlist(GS))
			# GET THE GENOTYPE CALLS OUT
			GD <- d[1:(length(d[,x])-2),x]
			# MAKE DATA CHARACTER
			GD <- as.character(GD)
			#make a copy
			GC <- GD
		########GC <- strsplit(GC,"_")
			GC <- c(rbind(gsub(".*_","",GD),gsub("_.*","",GD)))
		#	for (UZ in 1:length(GC)) {if(length(GC[[UZ]]) < 2) {GC[[UZ]][2] <- GC[[UZ]][1]}}
			
			# SPLIT THE GENEOTYPING CALLS INTO
			GD <- unlist(GC)
		#	GD <- unlist(strsplit(GD,"_"))
		#	GR <- unique(na.omit(GD))
			GR <- CD	
			#################################
			# now find the SNPs
			SE <- strsplit(CD,"[A-Z]|-") # equivalent to GR
			SN <- strsplit(CD,"[0-9]") # equivalent to GR

			# check for missing
			for (r in 1:length(SN))
			{
			SN[[r]] <- SN[[r]][nchar(SN[[r]]) >0]
			}

			for (r in 1:length(SE))
			{
			SE[[r]] <- SE[[r]][nchar(SE[[r]]) >0]
			}


			# check for substitutions in the insertions
			for (r in 1:length(SE))
				if(any(duplicated(SE[[r]]))) 
				{
				#find the common position
				SEpos <- SE[[r]][duplicated(SE[[r]])]
				# find the insert				
				SEI <- SN[[r]][duplicated(SE[[r]])]
				# find the match
				SEM <- SN[[r]][SE[[r]]  %in% SEpos & !grepl("I",SN[[r]])]
				# substitue the insert
				SN[[r]][SE[[r]]  %in% SEpos & !grepl("I",SN[[r]])] <- gsub("I",SEM,SEI)
				# now remove the insert 
				SN[[r]] <- SN[[r]][!duplicated(SE[[r]])]
				
				}			


				SLO <- 1:length(SNPpos)
				#check for insertions
				if (any(grep("I",SN)))
					{
#					SNPpos <- AC(SE[grep("I",SN)][[1]])
					SNPpos <- sort(AC(unique(unlist(SE[grep("I",SN)]))))
					if (DIR == "X16")
					{
					# FIND THE POSITIONS OF THE SNPS
					SLOC <-  AC(OFS)-SNPpos+1
					}
					# positive offset
					if (DIR == "X0")
					{
					# FIND THE POSITIONS OF THE SNPS
					SLOC <-  AC(OFS) -1 + SNPpos
					}
					}
			##################################
			# get the reference SNPs
			REF <- NULL
			for (t in 1:length(SLOC)){ REF <- c(REF,substr(GENseq,SNPpos[t],SNPpos[t]))}


			# make a copy and then fill in the values 
				SG <- SN
			# add the reference
			for (SR in 1:length(SN))
				{
				if(length(SN[[SR]]) < length(SNPpos))
					{
					print("short")
					for (WER in (SNPpos[!(SNPpos %in% SE[[SR]])]))
						{
						SG[[SR]] <- c(SG[[SR]],REF[SNPpos == WER])
						print(SE[[SR]])
						SE[[SR]] <- c(SE[[SR]],WER)
						print(SE[[SR]])
						}
#					for (WER in (length(SN[[SR]])+1):length(SNPpos))
#						{SG[[SR]][WER] <- REF[WER]}	
					}
				}
			# now sort the list(s) in order, of SNPpos
			for(SR in 1:length(SE))
				{
				SG[[SR]] <- SG[[SR]][order(AC(SE[[SR]][1:length(SNPpos)]))]
				}

			#now add the reference to the insert
			for (SR in 1:length(SG))
				{
				for(FG  in (grep("I",SG[[SR]])))
					{
					SG[[SR]][FG] <- gsub("I",REF[FG],SG[[SR]][FG])
					}

				}

			# make numerical
			SM <- matrix(unlist(SG),ncol=length(SG))
			SM2 <- SM

		# call the reference "0"
		ALT <- NULL
		for(UZ in 1:(length(SM[1,]))) {SM[SM[,UZ] == REF,UZ] <- "0" }

			for(UZ in 1:(length(SM[,1]))) 
				{
				LS <- unique(SM[UZ,SM[UZ,] != "0"])
				LN <- 1:length(LS)
				for(ZU in LN)	{SM[UZ,SM[UZ,] == LS[ZU]] <- LN[ZU]}
				print(paste(UZ,LS,LN))
				ALT <- c(ALT,paste(LS,collapse=","))
				}

	#	for (UZ in length(GD):1) 
	#		{
	#		if(is.na(GD[UZ])) {GD <- append(GD,NA,UZ)}
	#		}


		GGF <- NULL
		for( UZ in seq(1,length(GD),2))
			{
			OU <- NULL
			if(is.na(GD[UZ])) {OU <- rep("./.",length(SLOC))} else 
			{OU <- paste(SM[,GD[UZ] == GR],SM[,GD[(UZ+1)] == GR],sep="/")}
			GGF <- cbind(GGF,OU)
			
		}


		# now reorientate GR to CD (the genotype calling to the count calling)
		REO <- NULL
		LoCD <- 1:length(CD)
		for(UZ in LoCD) {REO <- c(REO,LoCD[GR[UZ] == CD])}

		SM <- matrix(SM[,REO],ncol=length(REO))

		CCF <- NULL
		for( UZ in seq(1,length(CRindividual),dim(SM)[2])) # genotype
			{
			GO <- NULL
			for (ZU in 1:length(SLOC)) # allele
				{	
				GCAL <- NULL
				for(GT in sort(unique(SM[ZU,]))) {GCAL <- c(GCAL,sum(AC(CRindividual[UZ:(UZ+((dim(SM)[2]-1)))][SM[ZU,] == GT]))) }
				GCAL <- paste(unlist(GCAL),collapse=",")
				GO <- c(GO,GCAL)
				}
			CCF <- cbind(CCF,GO)
			}

# chnage the reverse direction back to genome orientated nucleotides
		# negative offset
		if (DIR == "X16")
		{
		REF <- chartr("ATGC","TACG",REF)
		ALT <- chartr("ATGC","TACG",ALT)
		}

	###################################
	# check for "N"
	if (any(grep("N",GS)) || length(GS) < 2) {return(NA)} else
	{
	# SET UP A MATRIX FOR DATA OUPTUT
	VCF <- matrix(ncol=9,nrow=length(SLOC))
	# NOW POPULATE THE TABLE WITH REFERENCE DATA

	VCF[,1] <- CHR
	VCF[,2] <- SLOC
	VCF[,3] <- "."
	VCF[,4] <- REF
	VCF[,5] <- ALT
	VCF[,6] <- "ALT"
	VCF[,8] <- "."
	VCF[,9] <- "GT"
	VCF <- cbind(VCF,matrix(paste(GGF,CCF,sep=":"),ncol=dim(GGF)[2]))
	}
return(VCF)
}

OUT <- NULL

# single core
if (CORE == 1) {for (z in 1:(length(d)-1)){OUT <- rbind(OUT,makeVCF(z))}}


if (CORE >1) 
  {
  library(foreach)
  OUT <- NULL
  cl <- parallel::makeCluster(CORE)
  doParallel::registerDoParallel(cl)
  OUT <- foreach(i=(1:(length(d)-1)), .combine='rbind') %dopar% makeVCF(i)
}



OUTvcf <- as.data.frame(OUT)
colnames(OUTvcf) <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",as.character(d$Original[1:(length(d$Original)-2)]))

OUTvcf$FORMAT <- "GT:DPR"
OUTvcf$QUAL <- 100
OUTvcf[,2] <- AC(OUTvcf[,2])
#OUTvcf[,1] <- AC(OUTvcf[,1])

OUTvcf2 <- OUTvcf[order(OUTvcf[,1],OUTvcf[,2]),]
OUTvcf2 <- OUTvcf2[OUTvcf2$ALT != "",]
OUTvcf2 <- OUTvcf2[OUTvcf2$REF != "",]
OUTvcf2 <- OUTvcf2[!is.na(OUTvcf2$POS),]

colnames(OUTvcf2)[1] <- "##fileformat=VCFv\n#CHROM"
write.table(OUTvcf2,file=args[3],sep="\t",row.names=F,quote=F)




