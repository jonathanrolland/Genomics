# Script by D. Schluter and J. Rolland
# Convert vcf genotypes to snp matrices to be used in other programs
# using the file FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.vcf.gz"
# FirstFilter.GATK.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.vcf.gz contains all the chromosoms and scaffolds
# We remove the individual Ot021114-E15B5E that is YY from the vcf, named HI.4573.002.NS_Adaptor_23.Ot021114-E15B5E in the vcf and from the files with the populations
vcftools --remove-indv HI.4573.002.NS_Adaptor_23.Ot021114-E15B5E --gzvcf FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.vcf.gz --out FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.noYY.vcf --recode --recode-INFO-all

# compress
bgzip FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.noYY.vcf.recode.vcf

# ------------------------------------------------------------------------------------
# Generate a SNP matrix using snpStats

#FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.noYY.vcf.recode.vcf.gz
#GCA_018296145.1_Otsh_v2.0_genomic_wMito.gff

gzip -cd ~/Dropbox/1.en_cours/SticklebacksSalmons/121_chinook.v3/0.Kris_vcf_file_version3/FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.noYY.vcf.recode.vcf.gz | grep contig

# We extract the interesting part
# chr 1 -34 go from CM031199.1 to CM031232.1
# We copy-paste this information in scaffolds_chinook.txt
# create an index for the .vcf

tabix FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.noYY.vcf.recode.vcf.gz

#in R
library(snpStats)
library(data.table)
library(VariantAnnotation)
library(vcfR)
library(vroom)
library(data.table)
library(dplyr)

path<-"YOUR_PATH"
setwd(path)

#Now we create vcfs for each chromosoms, based on a small file built from the vcf file containing just scaffolds names.
scaff<-read.table("scaffolds_chinook.txt")$V1
for( i in 1:34 ){
system(paste0("bcftools view ",path,"FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.noYY.vcf.recode.vcf.gz --regions ", scaff[i]," | gzip > ",path,"chr",i,".chinook.var.vcf.gz"))
}

#However it should be run by chromosom to go faster
for( i in 1:34){
    print(i)
    vcf<-read.vcfR(file= paste0("chr",i,".chinook.var.vcf.gz") )
    snpmatrix <- extract.gt(vcf, element = 'GT')
    snpmatrix[which(snpmatrix == "1/1", arr.ind = TRUE)]<-2
 	snpmatrix[which(snpmatrix == "0/1", arr.ind = TRUE)]<-1
 	snpmatrix[which(snpmatrix == "1/0", arr.ind = TRUE)]<-1
 	snpmatrix[which(snpmatrix == "0/0", arr.ind = TRUE)]<-0
 	snpmatrix[which(is.na(snpmatrix), arr.ind = TRUE)]<-"NA"

    #keep only the SNP with QUAL>=20
    print(dim(snpmatrix))
    snpmatrix <- snpmatrix[which(as.numeric(vcf@fix[,6])>=20),]
    print(dim(snpmatrix))
    fwrite(cbind(rownames(snpmatrix),snpmatrix),file=paste0("chr",i,".snpmatrix_chinook.csv"),sep = ",")
}

##############################
# COMPUTE FREQUENCIES
##############################

# Random effects regressions of SNP allele frequencies on latitude
# River is the random effect

path<-"YOUR_PATH_INPUT"
path_output<-"YOUR_PATH_OUTPUT"
setwd(path)

#load latitude data
sampData <- read.csv("latitude_data.csv", sep=";")
scaff<- read.csv("scaffolds_chinook.txt")

for( i in 1:34){

snpmat<-data.frame(fread(file=paste0("chr",i,".snpmatrix_chinook.csv"),sep=","))
rownames(snpmat)<-snpmat[,1]
snpmat<-snpmat[,-1]

#transpose in order to have samples in the first column
snpmat<-t(snpmat)

fishNamesGeno <- sub(".*([O][t].*)", "\\1", rownames(snpmat))

# Check that samples are in the same order in the snpmat the sampData file
rownames(sampData)<-sampData$fishname
sampData<-sampData[fishNamesGeno,]
all(fishNamesGeno == sampData$fishname)

riverLatitude <- paste(sampData$river, sampData$latitude, sep = "|")

# Split genotype data frame by river-latitude
list_snpmat <- split(data.table(snpmat), riverLatitude)

# Total number of alleles by river summed across fish (multiplied the genotypes by 2 to get total alleles except for mtDNA)
totAlleles <- lapply(list_snpmat, function(x){ totAlleles <- apply(x, 2, function(x){2 * length(x[!is.na(x)])}) })

#paste all table together and write
totAllelesByRiver <- data.frame(do.call("rbind", totAlleles))
fwrite(totAllelesByRiver, file = paste0(path_output,"chr",i,".totAllelesByRiver.csv"), sep = ",", col.names = TRUE, row.names = TRUE)

# Sum of ALT alleles by river, summed across fish (divide by totAllelesByRiver to get allele frequencies)
nALTalleles <- lapply(list_snpmat, function(x){ nALTalleles <- apply(x, 2, sum, na.rm = TRUE)})

#paste all table together and write
nALTallelesByRiver <- data.frame(do.call("rbind", nALTalleles))
fwrite(nALTallelesByRiver, file = paste0(path_output,"chr",i,".nALTallelesByRiver.csv"), sep = ",", col.names = TRUE, row.names = TRUE)

# compute allele frequencies and write
freq <- nALTallelesByRiver/totAllelesByRiver

fwrite(freq, file = paste0(path_output,"chr",i,".alleleFreqByRiver.csv"), sep = ",", col.names = TRUE, row.names = TRUE)
}

##############################
# COMPUTE BGLMER
##############################

# Calculate allele frequencies at each site
# Carry out random effects linear regression of allele frequency by latitude

# in R
library(data.table)
library(blme)

path_output<-"PATH_OUTPUT"

for(i in 1:34){

# Read nREF and nALT
nalt <- as.data.frame(fread(paste0(path_output,"chr",i,".nALTallelesByRiver.csv")))
row.names(nalt) <- nalt[,1]
nalt <- nalt[,-1]
ntot <- as.data.frame(fread(paste0(path_output,"chr",i,".totAllelesByRiver.csv")))
row.names(ntot) <- ntot[,1]
ntot <- ntot[,-1]
nref <- ntot - nalt

latitude_init <- as.numeric(sub(".*[\\|](.*)", "\\1", rownames(nalt)))
river_init <- rownames(nalt)

# compute allele frequencies as proportions to remove the snp with MAF<0.05
# Sum REF and ALT alleles across populations
sumref <- apply(nref, 2, sum)
sumalt <- apply(nalt, 2, sum)
sumtot <- sumref + sumalt

# Drop snps having too low a minor allele frequency
# It looks like there are still a few in the data set

MAF.crit <- 0.05 # minor allele frequency threshold

freq <- sumalt / sumtot
#table(freq < MAF.crit | freq > 1 - MAF.crit) #for example: this deletes 212 out of 68981 snps for the CHR34

MAF <- sapply(freq, function(x){min(x, 1-x)})
nref <- nref[, MAF >= MAF.crit] 
nalt <- nalt[, MAF >= MAF.crit]

# The loop -- *** testing here only on the first 1000 snps ***
# give the table with the number of alt alleles and the number of ref alleles, individuals are in lines

npop<-dim(nalt)[1]

bglmer_table<-matrix(NA,length(colnames(nalt)),6)
bglmer_table[,1]<-colnames(nalt)
colnames(bglmer_table)<-c("snp","intercept","slope","pvalue","is.singular","fail.converge")

for (j in 1:dim(nalt)[2]){
	if(j%%100==0 | j==1)print(paste0(j,"/",dim(nalt)[2]))
	
	ref0 <- rep( rep(0, npop), nref[,j]) # we put a 0 for each individual with a ref allele in all populations
	alt1 <- rep( rep(1, npop), nalt[,j]) # we put a 1 for each individual with a alt allele in all populations
	river0 <- rep(river_init , nref[,j])
	river1 <- rep(river_init , nalt[,j])
	latitude0 <- rep(latitude_init , nref[,j])
	latitude1 <- rep(latitude_init , nalt[,j])
	p <- c(ref0, alt1)
	river <- c(river0, river1)
	latitude <- c(latitude0, latitude1)
	
	cleantable<-data.frame(p=as.numeric(p), latitude=as.numeric(latitude), river=as.character(river))
	
	#run bglmer and summarize results
 	z<-try(bglmer(cleantable$p ~ cleantable$lat + (1|cleantable$river), family = binomial(link = "logit")),silent=T)
 	if(!is.character(z)){
 	zSummary <- summary(z)
	z1 <- unname(coef(zSummary))
	fail.converge <- any(grepl("failed to converge", zSummary$optinfo$conv$lme4$messages))
	if(length(fail.converge) < 1) fail.converge <- FALSE
	bglmer_table[j,2:6]<-c(intercept = z1[1,1], slope = z1[2,1], P = z1[2,4], is.singular = isSingular(z), fail.converge = fail.converge) # glmer and bglmer version
 	}}

fwrite(bglmer_table, file = paste0(path_output,"bglmer.results.chr",i,".csv"), sep = ",", col.names = TRUE)
}



