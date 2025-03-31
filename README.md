# Genomics Chinook

This folder contains all the codes necessary to run the genomic analysis of the article "Genetic basis of latitudinal adaptation in Chinook salmon".

This includes the pipeline to go from the raw fastq files "XXX_R1.fastq.gz" and "XXX_R2.fastq.gz" available here: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA694998 to the final results of the analysis.

00.script.K.Christensen.pipeline.txt describes the pipeline to align and filter the reads and generate the vcf file FirstFilter.GATK.gt3.vcftools.biallele.mm0.9.maf0.05.meanDP5-200 used for further analyses. This script was run on cluster using SLURM.
Using the reference genome here: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/296/145/GCF_018296145.1_Otsh_v2.0/GCF_018296145.1_Otsh_v2.0_genomic.gff.gz

01.script.make.snp.matrices.and.run.bglmer.r explains how to make SNP matrices from the vcf and run the BGLMER analysis.

02.script.PCA.and.figures.r explains how to run PCA analyses and produce the figures.

03.script.admixture.analysis.r explains how to run the ADMIXTURE analyses.

04.script.BAYPASS.r describes how to run the BAYPASS analyses.

05.simplify.gff.r produces a simplified gff to merge information from the gff and the BAYPASS/BGLMER analyses.

06.script.compute.mean.results.per.gene.and.sliding.window.r compute the mean -log10pvalue and Bayes Factor per gene and per sliding window of 50 SNP.

07.script.analyses.results.BAYPASS.and.BGLMER.r script computing the results found in the main text testing whether each region fall in the top 5%.

08.Liftover analysis.merging.our.137indiv.with.160indiv.thompson.txt merge the two VCFs from our analysis and Thompson et al. 2020 analysis.

09.Comparison.ALASKA.and.CALIFORNIA.after.liftover.txt describes how to merge datasets and run analyses in order to analyse the spring-run alleles (using california samples from Thompson et al.) across the latitudinal gradient.

