# Script by J. Rolland 
###########
#Simplify the GFF

library(ape)
setwd("YOUR_PATH")
gff<-read.gff("GCA_018296145.1_Otsh_v2.0_genomic_wMito.gff") # annotation file from the reference genome

all_info_genes<-as.character(gff[which(gff[,3]=="gene"),9])
chr_gene<-as.character(gff[which(gff[,3]=="gene"),1])
start_gene<-as.numeric(gff[which(gff[,3]=="gene"),4])
end_gene<-as.numeric(gff[which(gff[,3]=="gene"),5])
strand_gene<-as.character(gff[which(gff[,3]=="gene"),7])

name_gene1<-rep(NA,length(all_info_genes))
for(k in 1:length(all_info_genes)) name_gene1[k]<-strsplit(strsplit(all_info_genes[k],split="Name=")[[1]][2],split=";")[[1]][1]

name_gene<-paste(name_gene1,chr_gene,start_gene,end_gene,sep="_")
simple_gff<-cbind(name_gene,chr_gene,start_gene,end_gene,strand_gene, all_info_genes)
colnames(simple_gff)<-c("gene","chr","start","end","strand","info")

#####
# Include here the regulatory regions
# rank by genes and by chromosom
# Give a name to the region  "reg_gene1_gene2" defined by the positions end+1 et start-1 of the following gene

regulatory_regions<-matrix(NA,dim(simple_gff)[1]+1,dim(simple_gff)[2])
colnames(regulatory_regions)<-c("gene","chr","start","end","strand","info")

for (i in 0:dim(regulatory_regions)[1]){
print(i)
if(i==0){ # for the very beginning of chr1
regulatory_regions[1,"gene"]<-paste(paste("REG",name_gene1[i+1],sep="|"),simple_gff[1,"chr"],1,as.numeric(simple_gff[i+1,"start"])-1,sep="_")
regulatory_regions[1,"chr"]<-simple_gff[1,"chr"]
regulatory_regions[1,"start"]<-1
regulatory_regions[1,"end"]<-as.numeric(simple_gff[i+1,"start"])-1
}
if(i==47189){# for the end of the last scaffold
regulatory_regions[i+1,"gene"]<-paste(paste("REG",name_gene1[i],sep="|"),simple_gff[i,"chr"],as.numeric(simple_gff[i,"end"])+1,1000000,sep="_")
regulatory_regions[i+1,"chr"]<-simple_gff[i,"chr"]
regulatory_regions[i+1,"start"]<-as.numeric(simple_gff[i,"end"])+1
regulatory_regions[i+1,"end"]<-1000000 #arbitrary high value
}
if(i > 0 & i < 47189) {
regulatory_regions[i+1,"gene"]<-paste(paste("REG",name_gene1[i],name_gene1[i+1],sep="|"),simple_gff[i,"chr"],as.numeric(simple_gff[i,"end"])+1,as.numeric(simple_gff[i+1,"start"])-1,sep="_")
regulatory_regions[i+1,"chr"]<-simple_gff[i,"chr"]
regulatory_regions[i+1,"start"]<-as.numeric(simple_gff[i,"end"])+1
regulatory_regions[i+1,"end"]<-as.numeric(simple_gff[i+1,"start"])-1
}
}#for

simple_gff_reg<-rbind(simple_gff,regulatory_regions)

fwrite(simple_gff_reg, file = "simple_gff.csv", sep = ",", col.names = TRUE, row.names = TRUE)
