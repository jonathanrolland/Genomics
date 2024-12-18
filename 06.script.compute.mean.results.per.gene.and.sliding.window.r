# Script by J. Rolland
# Merge the results of BAYPASS and the GFF and compute average results per gene and sliding window of 50 SNPs
# in order to do this I use a simplified version of the GFF computed in script 05.simplify.gff keeping only the gene (noted as gene in the gff), chr (chromosome or scaffold), start (position), end (position), strand and info.
#in R

library("data.table")
path_baypass<-"/Baypass_analysis/"

bf<-fread(file = paste0(path_baypass,"cache/baypass_env_model3_summary_betai_reg.out",h=T)$BF.dB.
snp<-fread(file = paste0(path_baypass,"cache_output_vcf2baypass_input_baypass/baypass_files3.snp",h=F)

#Compute the mean per gene and add the information of the gene on all positions in the genome
#for baypass

simple_gff<-read.table(file = "simple_gff.csv", sep = ",", h=TRUE)

chrnames<-unique(snp[,1])[1:34]
gene<-rep(NA,dim(snp)[1])

for(i in 1:dim(snp)[1]){
smallgff<-simple_gff[which(simple_gff$chr==snp[i,1]),]
   if(length(which(smallgff$end>=as.numeric(snp[i,2])))>0) {end<-which(smallgff$end>=as.numeric(snp[i,2])
   if(length(which(smallgff$start<=as.numeric(snp[i,2])))>0) {start<-which(smallgff$start<=as.numeric(snp[i,2]))
   gene[i]<-as.character(smallgff[intersect(end,start),"gene"])
  if(i%%5000==0) print(i)
   }}}
snp$gene<-gene
snp$bf<-bf

fwrite(snp,file=paste0(path_baypass,"snp_bf3.csv"),sep=",")

# Rank genes based on their BF
nb_genes<-length(unique(gene))
bf_per_gene<-rep(NA,nb_genes)
names(bf_per_gene)<-unique(gene)
for (i in 1:nb_genes){
    bf_per_gene[i]<-mean(snp$bf[which(snp$gene==names(bf_per_gene)[i])],na.rm=T)
    if(i%%5000==0) print(i)
}

genes_ranked_baypass<-bf_per_gene[order(bf_per_gene,decreasing=F)]
genes_ranked_baypass2<-cbind(names(genes_ranked_baypass),genes_ranked_baypass)
colnames(genes_ranked_baypass2)<-c("gene","BF")
rownames(genes_ranked_baypass2)<-NULL

fwrite(genes_ranked_baypass2,file=paste0(path_baypass,"genes_ranked_baypass3.txt"),sep=",",col.names=T)

####
# Same for BGLMER results 
# compute the average of the -log10(pvalue) and rank the genes

path_bglmer<-"YOUR_PATH"

list_files<-list.files(path_bglmer)
list_results<-list_files[grep(list_files,pattern="bglmer.results.chr")]

results<-fread(file = paste0(path_bglmer,list_results[1]), sep = ",", h=T)

for(i in 2:length(list_results)){
results<-rbind(results, fread(file = paste0(path_bglmer,list_results[i]), sep = ",", h=T))
}

#remove snps that have failed to converge
results<-results[-which(results[,6]==1),]

#compute chr name and position
snp<-cbind(sub("(.*)[\\_].*", "\\1", results$snp), as.numeric(sub(".*[\\_](.*)", "\\1", results$snp)))
colnames(snp)<-c("chr","position")

simple_gff<-read.table(file = "simple_gff.csv", sep = ",", h=TRUE)

chrnames<-unique(snp[,1])[1:34]
gene<-rep(NA,dim(snp)[1])

for(i in 1:dim(snp)[1]){
smallgff<-simple_gff[which(simple_gff$chr==snp[i,1]),]
   if(length(which(smallgff$end>=as.numeric(snp[i,2])))>0) {end<-which(smallgff$end>=as.numeric(snp[i,2]))
   if(length(which(smallgff$start<=as.numeric(snp[i,2])))>0) {start<-which(smallgff$start<=as.numeric(snp[i,2]))
   gene[i]<-as.character(smallgff[intersect(end,start),"gene"])
  if(i%%5000==0) paste(print(i),"/",dim(snp)[1])
   }}}
   
snp<-as.data.frame(snp)
snp$gene<-gene
snp$pvalue<-results$pvalue

fwrite(snp,file=paste0(path_bglmer,"snp_bglmer3.csv"),sep=",")


###
# Rank the bglmer results for each gene
snp<-fread(file=paste0(path_bglmer,"snp_bglmer3.csv"),sep=",")
gene<-snp$gene
nb_genes<-length(unique(gene))
pvalue_per_gene<-rep(NA,nb_genes)
names(pvalue_per_gene)<-unique(gene)
for (i in 1:nb_genes){
#mean of the -log10pvalue
    pvalue_per_gene[i]<-mean(-log10(snp$pvalue[which(snp$gene==names(pvalue_per_gene)[i])]),na.rm=T)
    if(i%%5000==0) print(i)
}

#order gene per ranked pvalues from smallest to largest
genes_ranked<-pvalue_per_gene[order(as.numeric(pvalue_per_gene),decreasing=T)]
genes_ranked2<-cbind(names(genes_ranked),genes_ranked)
colnames(genes_ranked2)<-c("gene","minuslog10pvalue")
rownames(genes_ranked2)<-NULL

fwrite(genes_ranked2,file=paste0(path_bglmer,"genes_ranked_bglmer3_log10.csv"),sep=",",col.names=T)


########
#Compute results for sliding windows

#the two results files are for all snps in the genome
path_baypass<-"/Baypass_analysis/"
path_bglmer<-"/1.allele_frequencies/"
path_figure<-"/figures/"


#for baypass 
library(ggmanh)
library(data.table)

baypass<-data.frame(fread(file=paste0(path_baypass,"snp_bf3.csv"),sep=",",h=T))
colnames(baypass)<-c("chr","position","gene","bf")

baypass_windows50<-matrix(NA,length(baypass$bf),2)
colnames(baypass_windows50)<-c("gene","meanpval_windows50")

for( i in 1:(length(baypass$bf)-24)){
baypass_windows50[i+24,2]<-mean(baypass$bf[i:(i+49)])
baypass_windows50[i+24,1]<-baypass$gene[i+24]
#print(i)
}

#order the table from high bf to low bf
baypass_windows50_ordered<-baypass_windows50[order(as.numeric(baypass_windows50[,2]),decreasing=T),]

#keep the best window per gene.
gene_bestwindow_baypass<-baypass_windows50_ordered[-which(duplicated(baypass_windows50_ordered[,1])),]

write.table(gene_bestwindow_baypass,file=paste0(path_baypass,"gene_bestwindow_baypass.txt"))

####
# for bglmer

bg<-data.frame(fread(file=paste0(path_bglmer,"snp_bglmer3.csv"),sep=",",h=T))

bg$pvalue <- -log10(as.numeric(bg$pvalue))

bglmer_windows50<-matrix(NA,length(bg$pvalue),2)
colnames(bglmer_windows50)<-c("gene","meanpval_windows50")

for( i in 1:(length(bg$pvalue)-24)){
bglmer_windows50[i+24,2]<-mean(as.numeric(bg$pvalue[i:(i+49)]))
bglmer_windows50[i+24,1]<-bg$gene[i+24]
#print(i)
}

#order the table from high log10pvalue to low log10pvalue
bglmer_windows50_ordered<-bglmer_windows50[order(as.numeric(bglmer_windows50[,2]),decreasing=T),]

#keep the best window per gene.
gene_bestwindow_bglmer<-bglmer_windows50_ordered[-which(duplicated(bglmer_windows50_ordered[,1])),]

write.table(gene_bestwindow_bglmer,file=paste0(path_bglmer,"gene_bestwindow_bglmer.txt"))

