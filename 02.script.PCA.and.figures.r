# From J. Rolland and Y. Assassa scripts
# Also followed the tutorial here: https://speciationgenomics.github.io/pca/

############ Plink/plink2
#the file FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.noYY.vcf.recode.vcf contains 8055508  variants and 137 individuals

#transform this large file in Plink format ped map bed and fam 
./plink   --vcf FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.noYY.vcf.recode.vcf  --recode --out Salmons3 --set-missing-var-ids @:# --double-id --allow-extra-chr

#8,055,508 variants and 137 people pass filters and QC.

#Prune the snp to create a file with only independent snps (according to the LD): this first line create a file with all the snp to remove
#From the description here: http://zzz.bwh.harvard.edu/plink/summary.shtml ; The parameters for --indep are: window size in SNPs (e.g. 50), the number of SNPs to shift the window at each step (e.g. 5), the VIF threshold. 
#The VIF is 1/(1-R^2) where R^2 is the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously. 
#That is, this considers the correlations between SNPs but also between linear combinations of SNPs. 
#A VIF of 10 is often taken to represent near collinearity problems in standard multiple regression analyses (i.e. implies R^2 of 0.9). 
#A VIF of 1 would imply that the SNP is completely independent of all other SNPs. 
#Practically, values between 1.5 and 2 should probably be used; particularly in small samples, if this threshold is too low and/or the window size is too large, too many SNPs may be removed.

./plink --file Salmons3 --indep 50 5 2 --set-missing-var-ids @:# --double-id --allow-extra-chr
  
#this removes 6299349 of 8055508 variants removed. -> this lead to only 1,756,159 remaining variants.
./plink --file Salmons3 --extract plink.prune.in --make-bed --out Salmons3_pruned --set-missing-var-ids @:# --double-id --allow-extra-chr

# Run PCA
./plink --bfile Salmons3_pruned --pca --allow-extra-chr --out Salmons3_pruned_PCA

#in R
library(nloptr)
library(car)
library(factoextra)
library(magrittr)
library(ade4)
library(factoextra)
library(magrittr)
library(ggplot2)
library(ggrepel)

setwd("YOUR_PATH")

##### Import eigen.vector file
PCA<- read.table("Salmons3_pruned_PCA.eigenvec", header = FALSE)[,-1]

#### Define Axis (=principal component) of the PCA
Axis1=PCA$V3
Axis2=PCA$V4
Axis3=PCA$V5

#### Add column of the percentage of variance explained by first Axis of the PCA
eig.val=read.table("Salmons3_pruned_PCA.eigenval", header=FALSE)
eig.val$variance.percent<-(eig.val$V1/sum(eig.val$V1)*100)

#not sure about this
#this gives the weight of each vector for the first 20 axes
#12.149227  8.907710  7.820811  6.986174  6.499243  5.795101  5.470360  5.221522  4.941647  4.707454  4.016393  3.641141  3.353966  3.248794  3.082978  2.909340  2.849329  2.823529  2.816320  2.758962

sampData <- read.table("latitude_data.csv", sep=';',h=T)

#load latitude data and remove individual Ot021114-E15B5E that is YY"
sampData <- sampData[-dim(sampData)[1],]
rownames(sampData)<-sampData$fishname

salmon_indiv<-read.table(file="Salmons3_pruned.fam",sep=" ")$V1
clean_names<-unlist(lapply(salmon_indiv,function(x){
split<-strsplit(x,split="\\.")[[1]]
split[length(split)]
}))

sampData<-sampData[clean_names,]
Latitude<-sampData$latitude
Longitude<-sampData$longitude
mid<-mean(Latitude)

sampData$state[which(sampData$state=="AL")]<-"AK"

#### Make the plot of the PCA for the regions

pdf(file= "plot PCA regions_latitude2.pdf",height=15, width=10)
ggplot(PCA,aes(x= V4, y= V3*-1, label=sampData$state))+geom_point(aes(color=Latitude),pch=16, cex=3)+
  scale_color_gradient2(midpoint = mid, low =  "#FC4E07" , mid = "#E7B800",high ="#00AFBB")+
  #geom_text_repel(min.segment.length = 0.5, seed = NA, box.padding = 0.5, max.overlaps = 8, aes(color=Latitude), force = 10,  max.time = 1, max.iter = Inf,)+
  geom_vline(xintercept = 0, linetype ="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="PC2 (8.91%)")+
  scale_y_continuous(name="PC1 (12.15%)")+
  ggtitle("Individuals - PCA")+
  theme(panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"))
dev.off()

  #### Make the plot of the PCA for the individuals
  pdf(file= "plot PCA individuals_latitude.pdf",height=20, width=40)
ggplot(PCA,aes(x= V3, y= V4, label=sampData$samples))+geom_point(aes(color=Latitude),pch=16, cex=3)+
  scale_color_gradient2(midpoint = mid, low =  "#FC4E07" , mid = "#E7B800",high ="#00AFBB")+
  geom_text_repel(min.segment.length = 0.5, seed = NA, box.padding = 0.5, max.overlaps = 8, aes(color=Latitude), force = 10,  max.time = 1, max.iter = Inf,)+
  geom_vline(xintercept = 0, linetype ="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="PC1 (12.15%)")+
  scale_y_continuous(name="PC2 (8.91%)")+
  ggtitle("Individuals - PCA")+
  theme(panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"))
       dev.off() 
        
          #### Make the plot of the PCA for the rivers
  pdf(file= "plot PCA river_latitude.pdf",height=15, width=30)
ggplot(PCA,aes(x= V3, y= V4, label=sampData$river))+geom_point(aes(color=Latitude),pch=16, cex=3)+
  scale_color_gradient2(midpoint = mid, low =  "#FC4E07" , mid = "#E7B800",high ="#00AFBB")+
  geom_text_repel(min.segment.length = 0.5, seed = NA, box.padding = 0.5, max.overlaps = 8, aes(color=Latitude), force = 10,  max.time = 1, max.iter = Inf,)+
  geom_vline(xintercept = 0, linetype ="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="PC1 (12.15%)")+
  scale_y_continuous(name="PC2 (8.91%)")+
  ggtitle("Individuals - PCA")+
  theme(panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"))
       dev.off() 
        
