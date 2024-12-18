# Code from J. Rolland
# I followed the pipeline proposed here
# https://github.com/Bio-protocol/Genome-Environment_Association_Analyses_with_BayPass
# download here: https://www1.montpellier.inra.fr/CBGP/software/baypass/
# untar the file
# install homebrew
# install the gfortran compiler
#run > make clean all FC=gfortran  #in the folder /Documents/baypass_2.4/sources
#type the following command to locate the path for launching
#export PATH=$PATH:"/Documents/baypass_2.4/sources"
# simply type > g_baypass to run the program

PATH=$PATH:"/Baypass_analysis/baypass_2.4/sources"

############
# On the R console

#prepare the latitude file for Baypass
setwd("/Baypass_analysis")

#modify the table to remove the YY individual
table_info<-read.csv("latitude_data.csv",h=T,sep=";")
table_info<-table_info[-138,]
table_info$river[59]<-"Salmon_BC"

write.table(table_info[,2:3],file="individuals_pop.tsv",row.names=F)

## I have then manually removed the quotes in the text files

lat<-table_info$latitude
names(lat)<-table_info$river

pop_lat<-matrix(NA,length(unique(table_info$river)),2)
pop_lat[,1]<-unique(table_info$river)
pop_lat[,2]<-lat[unique(table_info$river)]

write.table(pop_lat,file="pop_lat.tsv",row.names=F)
# I then manually removed the quotes in the text files
write.table(t(pop_lat),file="pop_lat_baypass.txt",row.names=F)

##################
#prepare the Genetic file for the Baypass analysis
#Use the perl script vcf2baypass.pl of https://github.com/Bio-protocol/Genome-Environment_Association_Analyses_with_BayPass

# BE CAREFUL !!! replace manually spaces by tab delimitation (\t) in the population.tsv file!!!

# in the Bash consol
#Run vcf2baypass.pl (to make it work, vcf should be compressed with the .gz format)
zcat < FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200.recode.noYY.vcf.recode.vcf.gz  | perl vcf2baypass.pl individuals_pop.tsv cache_output_vcf2baypass_input_baypass/baypass_files3

############################
#Before running baypass make sure that the covariate file (with latitude) has the same order than the ".pop" and the allele count files.

pop.baypass<-read.table("/Baypass_analysis/cache_output_vcf2baypass_input_baypass/baypass_files3.pop")$V1
pop.lat<-read.table("/Baypass_analysis/pop_lat_baypass_with_pop_names.txt",h=T)
new.pop.lat<-pop.lat[pop.baypass]

write.table(new.pop.lat,file="/Baypass_analysis/pop_lat_baypass_with_pop_names_ordered.txt",row.names=F)
names(new.pop.lat)<-NULL
write.table(new.pop.lat,file="/Baypass_analysis/pop_lat_baypass_final.txt",row.names=F)

# in the Bash consol
npop=$(wc -l ~/Baypass_analysis/baypass_files3.pop | cut -d " " -f1)
#npop is a variable containing 29 populations separated by spaces

#create a cache folder in /baypass_2.4/sources/
mkdir /baypass_2.4/sources/

cd /baypass_2.4/sources/

#launch the analysis of the core model to generate the covariance matrix # this takes about 15h
./g_baypass -npop $npop -gfile /Baypass_analysis/baypass_files3.txt -outprefix cache/baypass_core3 -nthreads 20

# I don't know if we should rescale the latitude, I don't think so, so I removed the -scalecov argument # the omega matrix is coming from the core model output (previous command line)
./g_baypass -npop $npop -gfile /Baypass_analysis/baypass_files3.txt -efile ~/Baypass_analysis/pop_lat_baypass_final.txt  -omegafile cache/baypass_core3_mat_omega.out -outprefix cache/baypass_env_model3 -nthreads 20

## Plot population correlation from BayPass covariance matrix
# in R
library(tidyverse)
library(corrplot)
select=dplyr::select

setwd("/Baypass_analysis/cache")

my.cov <- as.matrix(read.table("baypass_core3_mat_omega.out"))
my.pop <- read.table("../cache_output_vcf2baypass_input_baypass/baypass_files3.pop")
rownames(my.cov)<-my.pop$V1
colnames(my.cov)<-my.pop$V1

#order the matrix by latitude for the plot
rt<-as.numeric(read.table("/Baypass_analysis/pop_lat_baypass_final.txt"))
my.cov2<-my.cov[order(rt),order(rt)]

my.corr <- cov2cor(my.cov2)
png("core_matrix3_omega.png", 750, 750)
corrplot(my.corr, mar=c(2,1,2,2)+0.1, main=expression("Correlation map based on"~hat(Omega)))
dev.off()




























