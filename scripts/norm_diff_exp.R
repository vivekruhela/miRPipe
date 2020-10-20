#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))

data_conditions <- read.csv("data/sample_list.csv")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	 stop("Please pass all the arguments in the following order: norm_diff_exp.R <rna_type>", call.=FALSE)
}

rna_type <- args[1]
if (rna_type == 'miRNA'){
	SEQ_DIR <- file.path(getwd(),'data')
	cat('Setting WORKDIR to:', SEQ_DIR, '\n' )
	setwd(SEQ_DIR)
	cts<- read.table('counts_final_1.csv', sep=",", header = T,stringsAsFactors=F)
} else {
	SEQ_DIR <- file.path(getwd(),'data','piRNA') 
	cat('Setting WORKDIR to:', SEQ_DIR, '\n' )
	setwd(SEQ_DIR)
	cts<- read.table('piRNA_raw_counts.csv', sep=",", header = T,stringsAsFactors=F)
}


cts1 = as.matrix(cts[,2:ncol(cts)])
gnames = as.matrix(cts[,1])
rownames(cts1) = gnames

samples = colnames(cts1)

# Experiment Design 
expDesign <- data.frame(row.names = colnames(cts1), condition = data_conditions$condition)


# Making DESeq2 object
dds <- DESeqDataSetFromMatrix(countData =  cts[,2:ncol(cts)], colData = expDesign, design = ~condition)
featureData <- data.frame(gene=rownames(cts1))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# Performing DESeq2 analysis
dds <- DESeq(dds)
normalized = counts(dds, normalized=TRUE)
res <- results(dds,contrast=c("condition","treated","untreated"))
#colData(dds)
joint<- cbind(cts[,1], res)

# Saving results
if (rna_type == 'miRNA') {
	write.csv(file='miRNA_DESeq_expression_data_controlIsUnreated.csv', x=joint)
	write.csv(file='miRNA_DESeq_norm_size_factor_data.csv', x=normalized)
	cat('Saving results in miRNADESeq_expression_data_controlIsUnreated.csv. ', '\n')

	# Identification of significantly expressed mi/pi/sno-RNA
	cat('DESeq2 analysis is complete. Now finding significantly expressed entries.....','\n')
	input_mir <- read.csv('miRNA_DESeq_expression_data_controlIsUnreated.csv')

	# Selecting differentially expressed miRNAs with p_adj <= 0.05
	output_mir <- (subset(input_mir,input_mir$padj <= 0.05,select=cts...1.))

	# Saving results
	write.table(file='miRNA_significantly_DE_mir.csv', x=output_mir,quote=FALSE,row.names=FALSE,col.names=FALSE)
	cat('Final the results are saved in miRNA_significantly_DE_mir.csv ')

} else {
	write.csv(file='piRNA_DESeq_expression_data_controlIsUnreated.csv', x=joint)
	write.csv(file='piRNA_DESeq_norm_size_factor_data.csv', x=normalized)
	cat('Saving results in piRNADESeq_expression_data_controlIsUnreated.csv. ', '\n')

	# Identification of significantly expressed mi/pi/sno-RNA
	cat('DESeq2 analysis is complete. Now finding significantly expressed entries.....','\n')
	input_mir <- read.csv('piRNA_DESeq_expression_data_controlIsUnreated.csv')

	# Selecting differentially expressed miRNAs with p_adj <= 0.05
	output_mir <- (subset(input_mir,input_mir$padj <= 0.05,select=cts...1.))

	# Saving results
	write.table(file='piRNA_significantly_DE_mir.csv', x=output_mir,quote=FALSE,row.names=FALSE,col.names=FALSE)
	cat('Final the results are saved in piRNA_significantly_DE_mir.csv ')

}
	

