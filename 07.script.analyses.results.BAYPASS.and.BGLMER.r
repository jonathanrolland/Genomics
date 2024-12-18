
####### MAIN RESULTS
# See if candidates genes are ranking in the top 5% of the regions (mean pvalue or sliding window 50snp)
# With mean -log10() pvalues from bglmer 
# and mean Bayes factor (BF) from Baypass

library("data.table")
path_baypass<-"/Baypass_analysis/"
path_bglmer<-"/1.allele_frequencies/"

baypass<-read.table(file=paste0(path_baypass,"genes_ranked_baypass3.txt"),sep=" ",h=T)
bg<-data.frame(fread(file=paste0(path_bglmer,"genes_ranked_bglmer3_log10.csv"),sep=",",h=T))

baypass<-baypass[,-1]

bg<-bg[order(bg[,"minuslog10pvalue"],decreasing=T),]
baypass<-baypass[order(baypass[,"BF"],decreasing=T),]

gene_bestwindow_baypass<-read.table(file=paste0(path_baypass,"gene_bestwindow_baypass.txt"),header=T)
gene_bestwindow_bglmer<-read.table(file=paste0(path_bglmer,"gene_bestwindow_bglmer.txt"),header=T)

#"clock1b"= clocka,"clock1a"=LOC112251141, "rock1","vgll3","six6"=LOC112233349
candidate_genes<-c("clocka","LOC112251141","rock1","vgll3","LOC112233349")

ranking_summary_list<-list()

for(i in 1:length(candidate_genes)){
    print("########################")
    print(candidate_genes[i])
    genes_bg<- grep(bg[,1], ignore.case =T,pattern=candidate_genes[i])
    genes_baypass<-grep(baypass[,1], ignore.case =T, pattern=candidate_genes[i],)
    
    max_genes<-max(length(genes_bg), length(genes_baypass))
    cond <-which.max(c(length(genes_bg), length(genes_baypass)))
    if(cond==1)names_genes<-bg[genes_bg,1]
    if(cond==2)names_genes<-baypass[genes_baypass,1]
    #if(sum(cond ==2 | cond==c(1,2))>0)names_genes<-baypass[genes_baypass,1]
   
    ranking_summary_table<-matrix(NA,max_genes,4)
    rownames(ranking_summary_table)<-names_genes
    colnames(ranking_summary_table)<-c("bglmer_mean-log10pvalue(%)","bglmer_50SNPwindow(%)","baypass_mean-log10pvalue(%)","baypass_50SNPwindow(%)")
    
    if( max_genes>0){for (j in 1:max_genes){
        
        if(sum(bg[,1]==names_genes[j],na.rm=T)>0)  	ranking_summary_table[j,1] 	<- 	round((which(bg[,1]==names_genes[j])/dim(bg)[1])*100,2)
        if(sum(bg[,1]==names_genes[j],na.rm=T)>0) 	ranking_summary_table[j,2]	<- 	round((which(gene_bestwindow_bglmer[,1]==names_genes[j])/dim(gene_bestwindow_bglmer)[1])*100,2)
    
        if(sum(baypass[,1]==names_genes[j],na.rm=T)>0) 	ranking_summary_table[j,3] 	<- 	round((which(baypass[,1]==names_genes[j])/dim(baypass)[1])*100,2)
        if(sum(baypass[,1]==names_genes[j],na.rm=T)>0) 	ranking_summary_table[j,4]	<-	round((which(gene_bestwindow_baypass[,1]==names_genes[j])/dim(gene_bestwindow_baypass)[1])*100,2)
    }}
    print(ranking_summary_table)
    ranking_summary_list[[i]]<-ranking_summary_table
}

##############################
# Plot interesting significant genes 
##############################

#The threshold quantiles at 95% of all snp are 3.249211 for the -log10(pvalue) of the bglmer and 2.914742 for the bayes factor of the baypass

#Example with clock1b
candidate_genes<-c("REG|tmem165|clocka_CM031199.1_30337921_30346550")

#Rank the regions according to their best SNP
listgene<-unique(snp$gene)
listmax<-rep(NA,length(listgene))
for( i in 1:length(listgene)){listmax[i]<-max(snp$bf[which(snp$gene==listgene[i])])}

for(analysis in c("baypass","bglmer") ){

if(analysis=="baypass")snp<-fread(file=paste0(path_baypass,"snp_bf3.csv"),sep=",",h=T)
if(analysis=="bglmer")snp<-fread(file=paste0(path_bglmer,"snp_bglmer3.csv"),h=T)

if(analysis=="baypass")colnames(snp)<-c("chr","position","gene","bf")
if(analysis=="bglmer")colnames(snp)<-c("chr","position","gene","pvalue")

snp<-as.data.frame(snp)
#if(analysis=="baypass")snp<-snp[-which(is.na(snp$bf)),]
if(analysis=="bglmer")snp<-snp[-which(is.na(snp$pvalue)),]

for(i in 1:length(candidate_genes)){
print(candidate_genes[i])
positions_gene<-which(snp$gene==candidate_genes[i])

#cut the file corresponding to the significant gene
snp_gene<-snp[positions_gene,]
chr<-snp_gene[1,1]

snp_chr<-snp[which(snp$chr==chr),]

min<-min(snp_gene$position) #identify the gene before and the gene after
max<-max(snp_gene$position) #which position is the closest but inferior

other_snp<-snp[-positions_gene,] # remove the positions inside the gene on the whole
other_snp_same_chr<-other_snp[which(other_snp$chr==chr),]

gene_before<-other_snp_same_chr[which.min(abs(other_snp_same_chr$position - min)),3]# find the gene before
gene_after<-other_snp_same_chr[which.min(abs(other_snp_same_chr$position - max)),3]# find the gene after

snp_print<-snp[c(which(snp$gene==gene_before), which(snp$gene==candidate_genes[i]),which(snp$gene==gene_after)),]

col<-rep("#999999",dim(snp_print)[1])

if(analysis=="baypass")significant_snp<-quantile(snp$bf,0.95,na.rm=T)
if(analysis=="baypass")col[which(snp_print$bf>significant_snp)]<-"#be2517"

if(analysis=="bglmer")significant_snp<-quantile(-log10(as.numeric(snp$pvalue)),0.95,na.rm=T)
if(analysis=="bglmer")col[which(-log10(as.numeric(snp_print$pvalue))>significant_snp)]<-"#be2517"

name_candidate_genes <-strsplit(candidate_genes[i],split="_")[[1]][1]
name_gene_before <- strsplit(gene_before,split="_")[[1]][1]
name_gene_after <-strsplit(gene_after,split="_")[[1]][1]

realname_chr<-as.numeric(sub("CM","",chr))-31198.1

if(sum(grep(name_candidate_genes,pattern="[|]"))>0) name_candidate_genes <- "Regulatory region"
if(sum(grep(name_gene_before,pattern="[|]"))>0)name_gene_before <- "Regulatory region"
if(sum(grep(name_gene_after,pattern="[|]"))>0)name_gene_after <- "Regulatory region"

par(mar=c(6,5,4,2))
if(analysis=="baypass")pdf(paste0("~/Dropbox/1.en_cours/SticklebacksSalmons/121_chinook.v3/figures/",candidate_genes[i],"_baypass.pdf"))
if(analysis=="baypass")plot(as.numeric(snp_print$position)*10^-6, as.numeric(snp_print$bf),pch=20,cex=1, col=col,cex.lab=1, cex.axis=0.75, cex.main=1.5,ylab="", xlab="",xlim=c(min(snp_print$position)*10^-6 - 0.003,max(snp_print$position)*10^-6 + 0.003), ylim=c(min(snp_print$bf),max(snp_print$bf))) 

if(analysis=="bglmer")pdf(paste0("~/Dropbox/1.en_cours/SticklebacksSalmons/121_chinook.v3/figures/",candidate_genes[i],"_bglmer.pdf"))
if(analysis=="bglmer")plot(as.numeric(snp_print$position)*10^-6, -log10(as.numeric(snp_print$pvalue)),pch=20,cex=1, col=col,cex.lab=1, cex.axis=0.75, cex.main=1.5,ylab="", xlab="",xlim=c(min(snp_print$position)*10^-6 - 0.003,max(snp_print$position)*10^-6 + 0.003), ylim=c(min(-log10(as.numeric(snp_print$pvalue))),max(-log10(as.numeric(snp_print$pvalue))))) 

if(analysis=="baypass")print(mean(as.numeric(snp_print$bf)))
if(analysis=="bglmer")print(mean(-log10(as.numeric(snp_print$pvalue))))
print(significant_snp)

#par(new=T)
#plot(as.numeric(red_dots$position)*10^-6, as.numeric(red_dots$bf),pch=20,cex=1, col="#be2517",cex.lab=1, cex.axis=0.75, cex.main=1.5,ylab="", xlab="",xlim=c(min*10^-6,max*10^-6),ylim=c(min(snp_print$bf),max(snp_print$bf))) 

#title(main = "Genomic location of variation associated with latitude within \n 249 kb surrounding the GREB1L/ROCK1 region in Chinook salmon")
title(xlab = paste0("Position on the chromosom ",realname_chr," (Mb)"), line = 4)

if(analysis=="baypass")title(ylab = "Bayes Factor", line = 3)
if(analysis=="bglmer")title(ylab = "-log10(P-value)", line = 3)

abline(v=min(snp_print$position)*10^-6, col="#d17805", lwd=1, lty=2)

abline(v=as.numeric(min)*10^-6, col="#be2517", lwd=1, lty=2)
abline(v=as.numeric(max)*10^-6, col="#be2517", lwd=1, lty=2)
abline(v=max(snp_print$position)*10^-6, col="#d17805", lwd=1, lty=2)

mtext(name_gene_before, side = 1, line = 2.9, at=(min(snp_print$position)*10^-6 + min*10^-6)/2, col="#d17805", cex=1)
mtext(name_candidate_genes, side = 1, line = 2.2, at=(min*10^-6 + max*10^-6)/2, col="#be2517", cex=1)
mtext(name_gene_after, side = 1, line = 2.9, at=(max*10^-6 + max(snp_print$position)*10^-6)/2, col="#d17805", cex=1)

dev.off()
print("number of significant snp, number of total snp, percent of significant snp")
if(analysis=="baypass")print(paste(length(which(snp_print$bf>significant_snp)),length(snp_print$bf),round((length(which(snp_print$bf>significant_snp))/length(snp_print$bf))*100,1)))
if(analysis=="bglmer")print(paste(length(which(-log10(as.numeric(snp_print$pvalue))>significant_snp)),length(snp_print$pvalue),round((length(which(-log10(as.numeric(snp_print$pvalue))>significant_snp))/length(snp_print$pvalue))*100,1)))

}
}

#######################################
# Plot the frequencies

#list clocka,rock1,vgll3 and six6 and their neighboring regions
list_gene<- c("clock1a_CM031203.1_65574213_65642249","clock1b_CM031199.1_30337921_30419389","greb1lrock1_CM031226.1_13415735_13619775","vgll3_CM031201.1_57876958_57939896","six6_CM031209.1_30434666_30489762") #I have changed the positions of the regions to include neighboring regions

pdf(file = paste0("~/Dropbox/1.en_cours/SticklebacksSalmons/121_chinook.v3/candidate_gene_baypass",list_gene[gene],".pdf"), width=10, height=15)
par(mfrow=c(5,2))
for (gene in c(1:5)){

gene_for_plot<- list_gene[gene]

results_bf <- data.frame(fread(paste0("~/Dropbox/1.en_cours/SticklebacksSalmons/121_chinook.v3/Baypass_analysis/all_chr_baypass.csv"), sep = ","))
 
if(gene==1)matrix_gene<-results_bf[c(which(results_bf$gene == "REG|LOC112251137|LOC112251141_CM031203.1_65574213_65581550"),which(results_bf$gene =="LOC112251141_CM031203.1_65581551_65632197"),which(results_bf$gene =="REG|LOC112251141|LOC112251822_CM031203.1_65632198_65642249")),]
if(gene==2)matrix_gene<-results_bf[c(which(results_bf$gene == "REG|tmem165|clocka_CM031199.1_30337921_30346550"),which(results_bf$gene =="clocka_CM031199.1_30346551_30407903"),which(results_bf$gene =="REG|clocka|nmu_CM031199.1_30407904_30419389")),]
if(gene==3)matrix_gene<-results_bf[c(which(results_bf$gene == "LOC112226774_CM031226.1_13415735_13496978"),which(results_bf$gene =="REG|LOC112226774|rock1_CM031226.1_13496979_13548901"),which(results_bf$gene =="rock1_CM031226.1_13548902_13619775")),]
if(gene==4)matrix_gene<-results_bf[c(which(results_bf$gene == "REG|LOC112239858|vgll3_CM031201.1_57876958_57887978"),which(results_bf$gene =="vgll3_CM031201.1_57887979_57892088"),which(results_bf$gene =="REG|vgll3|akap11_CM031201.1_57892089_57939896")),]
if(gene==5)matrix_gene<-results_bf[c(which(results_bf$gene == "REG|LOC112233348|LOC112233349_CM031209.1_30434666_30454759"),which(results_bf$gene =="LOC112233349_CM031209.1_30454760_30459932"),which(results_bf$gene =="REG|LOC112233349|LOC112233351_CM031209.1_30459933_30489762")),]


threshold95baypass<-2.914742

col<-rep("#0000000D",length(matrix_gene$bf))
col[which(as.numeric(matrix_gene$bf)>threshold95baypass)]<-"#BE251780"
names(col)<-matrix_gene$position

plot( matrix_gene$position, matrix_gene$bf, ylab="Bayes factor",xlab="Position",main= gene_for_plot, col=col, pch=16)

if(gene==1)abline(v=65581551,lty=2);abline(v=65632197,lty=2)
if(gene==2)abline(v=30346551,lty=2);abline(v=30407903,lty=2)
if(gene==3)abline(v=13496979,lty=2);abline(v=13548901,lty=2)
if(gene==4)abline(v=57887979,lty=2);abline(v=57892088,lty=2)
if(gene==5)abline(v=30454760,lty=2);abline(v=30459932,lty=2)

hit_set<-gene_for_plot

for (i in 1:length(hit_set)){

chrname<-strsplit(hit_set[i],split="_")[[1]][2]
chrnumber<-as.numeric(sub(x=chrname,pattern="CM",replacement=""))-31198.1
position_start<-as.numeric(strsplit(hit_set[i],split="_")[[1]][3])
position_end<-as.numeric(strsplit(hit_set[i],split="_")[[1]][4])
freq <- data.frame(fread(paste0("chr",chrnumber,".alleleFreqByRiver.csv"), sep = ","))
lat <- as.numeric(sub(".*[\\|](.*)", "\\1", freq[,1]))
position_freq<-sapply(colnames(freq),function(x)as.numeric(strsplit(x,split="_")[[1]][2]))

#select the columns in freq that correspond to the gene
col_selected<-intersect(which(position_freq>position_start),(which(position_freq<position_end)))

for(p in 1:length(col_selected)){
print(p)

if(length(which(is.na(freq[,col_selected[p]])))>0){
freq_plot <- freq[,col_selected[p]][-which(is.na(freq[,col_selected[p]]))]
lat_plot <- lat[-which(is.na(freq[,col_selected[p]]))]}else{
freq_plot <-  freq[,col_selected[p]]
lat_plot <- lat
}

position_freq_selected<-position_freq[col_selected][p]

#reorder all frequencies, with high latitude being the reference (frequency increasing to high latitude)
if( freq_plot[which(lat_plot == min(lat_plot))] > freq_plot[which(lat_plot == max(lat_plot))]) freq_plot<-sapply(freq_plot,function(x)1-x)

ss<-smooth.spline(lat_plot,freq_plot,df=5)
#if(P_position>=0.05){ 
plot(ss$x,ss$y,type="l",xlim=c(37,62), ylim=c(0,1),xlab="Latitude(°)",ylab="Frequency",main=(hit_set[i]),col=col[as.character(position_freq_selected)])

par(new=T)
if(p==length(col_selected))plot(ss$x,ss$y,type="l",xlim=c(37,62), ylim=c(0,1),xlab="Latitude(°)",ylab="Frequency",main=(hit_set[i]),col=col[as.character(position_freq_selected)]) # to end the loop

}
}
}
dev.off()


##########################
# Same but for bglmer

#I have manually changed the positions of the regions to include neighboring regions for the plot
list_gene<- c("clock1a_CM031203.1_65574213_65642249","clock1b_CM031199.1_30337921_30419389","greb1lrock1_CM031226.1_13415735_13619775","vgll3_CM031201.1_57876958_57939896","six6_CM031209.1_30434666_30489762") 

pdf(file = paste0("candidate_gene_bglmer_",list_gene[gene],".pdf"), width=10, height=15)
par(mfrow=c(5,2))
for (gene in c(1:5)){

gene_for_plot<- list_gene[gene]

results_bf <- data.frame(fread(paste0("all_chr_bglmer.csv"), sep = ","))
 
if(gene==1)matrix_gene<-results_bf[c(which(results_bf$gene == "REG|LOC112251137|LOC112251141_CM031203.1_65574213_65581550"),which(results_bf$gene =="LOC112251141_CM031203.1_65581551_65632197"),which(results_bf$gene =="REG|LOC112251141|LOC112251822_CM031203.1_65632198_65642249")),]
if(gene==2)matrix_gene<-results_bf[c(which(results_bf$gene == "REG|tmem165|clocka_CM031199.1_30337921_30346550"),which(results_bf$gene =="clocka_CM031199.1_30346551_30407903"),which(results_bf$gene =="REG|clocka|nmu_CM031199.1_30407904_30419389")),]
if(gene==3)matrix_gene<-results_bf[c(which(results_bf$gene == "LOC112226774_CM031226.1_13415735_13496978"),which(results_bf$gene =="REG|LOC112226774|rock1_CM031226.1_13496979_13548901"),which(results_bf$gene =="rock1_CM031226.1_13548902_13619775")),]
if(gene==4)matrix_gene<-results_bf[c(which(results_bf$gene == "REG|LOC112239858|vgll3_CM031201.1_57876958_57887978"),which(results_bf$gene =="vgll3_CM031201.1_57887979_57892088"),which(results_bf$gene =="REG|vgll3|akap11_CM031201.1_57892089_57939896")),]
if(gene==5)matrix_gene<-results_bf[c(which(results_bf$gene == "REG|LOC112233348|LOC112233349_CM031209.1_30434666_30454759"),which(results_bf$gene =="LOC112233349_CM031209.1_30454760_30459932"),which(results_bf$gene =="REG|LOC112233349|LOC112233351_CM031209.1_30459933_30489762")),]

matrix_gene$log10pvalue<- -log10(as.numeric(matrix_gene$pvalue))

threshold95bglmer<-3.249211 #obtained from quantile(-log10(as.numeric(snp$pvalue)),0.95,na.rm=T) in another script

col<-rep("#0000000D",length(matrix_gene$log10pvalue))
col[which(as.numeric(matrix_gene$log10pvalue)>threshold95bglmer)]<-"#BE251780"
names(col)<-matrix_gene$position

plot( matrix_gene$position, matrix_gene$log10pvalue, ylab="-log10pvalue",xlab="Position",main= gene_for_plot, col=col, pch=16) # pch = shapes, 

if(gene==1)abline(v=65581551,lty=2);abline(v=65632197,lty=2)
if(gene==2)abline(v=30346551,lty=2);abline(v=30407903,lty=2)
if(gene==3)abline(v=13496979,lty=2);abline(v=13548901,lty=2)
if(gene==4)abline(v=57887979,lty=2);abline(v=57892088,lty=2)
if(gene==5)abline(v=30454760,lty=2);abline(v=30459932,lty=2)

hit_set<-gene_for_plot

for (i in 1:length(hit_set)){

chrname<-strsplit(hit_set[i],split="_")[[1]][2]
chrnumber<-as.numeric(sub(x=chrname,pattern="CM",replacement=""))-31198.1
position_start<-as.numeric(strsplit(hit_set[i],split="_")[[1]][3])
position_end<-as.numeric(strsplit(hit_set[i],split="_")[[1]][4])
freq <- data.frame(fread(paste0("chr",chrnumber,".alleleFreqByRiver.csv"), sep = ","))
lat <- as.numeric(sub(".*[\\|](.*)", "\\1", freq[,1]))
position_freq<-sapply(colnames(freq),function(x)as.numeric(strsplit(x,split="_")[[1]][2]))

#select the columns in freq that correspond to the gene
col_selected<-intersect(which(position_freq>position_start),(which(position_freq<position_end)))

for(p in 1:length(col_selected)){
print(p)

if(length(which(is.na(freq[,col_selected[p]])))>0){
freq_plot <- freq[,col_selected[p]][-which(is.na(freq[,col_selected[p]]))]
lat_plot <- lat[-which(is.na(freq[,col_selected[p]]))]}else{
freq_plot <-  freq[,col_selected[p]]
lat_plot <- lat
}

position_freq_selected<-position_freq[col_selected][p]

if( freq_plot[which(lat_plot == min(lat_plot))] > freq_plot[which(lat_plot == max(lat_plot))]) freq_plot<-sapply(freq_plot,function(x)1-x)
ss<-smooth.spline(lat_plot,freq_plot,df=5)
plot(ss$x,ss$y,type="l",xlim=c(37,62), ylim=c(0,1),xlab="Latitude(°)",ylab="Frequency",main=(hit_set[i]),col=col[as.character(position_freq_selected)])
par(new=T)
if(p==length(col_selected))plot(ss$x,ss$y,type="l",xlim=c(37,62), ylim=c(0,1),xlab="Latitude(°)",ylab="Frequency",main=(hit_set[i]),col=col[as.character(position_freq_selected)]) # to end the loop
}
}
}
dev.off()
