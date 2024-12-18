############
# ADMIXTURE
############

# From Salmons3_pruned.bed file (file contained the unlinked sites), delete name of all chromosomes
awk '{$1=0;print $0}' ../PCA_salmon/Salmons3_pruned.bim > Salmons3_pruned2.bim

#Copy and paste the files .bed and .fam from the PCA_analysis

# Run ADMIXTURE for K (from 1 to 10)
for K in 1 2 3 4 5 6 7 8 9 10; \
do ./admixture --cv Salmons3_pruned2.bed $K |\
tee Salmons3_pruned.${K}.out;\
done

# Display the cross-validation values of each file which ends with .out
cat Salmons3*out | grep CV

##### Barplot of admixture results K=2

#in R
sampData <- read.table("latitude_data",h=T,sep=";")
sampData <- sampData[-dim(sampData)[1],]
rownames(sampData)<- sampData[,1]

#load coef of ancestry
tbl<-read.table("Salmons3_pruned2.2.Q")

#when running the analysis of admixture the individuals were in this order:
rownames(tbl)<- read.table("Salmons3_pruned2.fam",h=F)$V1
#write.table(tbl,"tbl.txt")

#sampData and names_order_vcf have the same size but we have to modify the order of sampData because the names are not matching completely
ordering<-matrix(NA,length(rownames(tbl)),2)
colnames(ordering)<-c("position_sampData","position_tbl")

for ( i in 1:length(rownames(tbl))){
ordering[i,1]<-i
ordering[i,2]<-grep(rownames(tbl),pattern=sampData[i,1])
}

tbl2<-as.matrix(tbl[ordering[,2],])

pdf(file="Plot ancestry_admixture (K=2).pdf",height=8,width=25)
par(mar=c(25,4,3,0))
barplot(t(tbl2), col=c("#FC4E07","#00AFBB"), space = c(0.5,0.5),
        names.arg = paste(sampData$samples,sampData$river,sampData$state, sep="_"), las = 2,
        ylab="Ancestry", main="Admixture software partitions variation in complete genome into two clusters", border=NA)
        
dev.off()
