# For MiRPipe expression data
SEQ_DIR <- file.path(getwd(),'data')

# For moRNA expression data
# SEQ_DIR <- file.path(getwd(),'data','moRNA','hairpin_counts') 


# For piRNA expression data
# SEQ_DIR <- file.path(getwd(),'data','piRNA','pirna_counts') 


# For snoRNA expression data
# SEQ_DIR <- file.path(getwd(),'data','snoRNA','snorna_counts') 

cat('Setting WORKDIR to:', SEQ_DIR, '\n' )
setwd(SEQ_DIR)

suppressPackageStartupMessages(library(DESeq2))

# Loading the MiRPipe merged gene count data
cts<- read.table('count_merged_1.csv', sep=",", header = T,stringsAsFactors=F)

# Loading the moRNA merged gene count data
#cts<- read.table('morna_count.csv', sep=",", header = T,stringsAsFactors=F)

# Loading the piRNA merged gene count data
#cts<- read.table('piRNA_raw_counts.csv', sep=",", header = T,stringsAsFactors=F)

# Loading the snoRNA merged gene count data
#cts<- read.table('snoRNA_raw_counts.csv', sep=",", header = T,stringsAsFactors=F)


cts1 = as.matrix(cts[,2:ncol(cts)])
gnames = as.matrix(cts[,1])
rownames(cts1) = gnames

samples = colnames(cts1)

# Experiment Design 
expDesign <- data.frame(row.names = colnames(cts1), condition = c(rep("treated",3),"untreated","untreated"))


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
write.csv(file='DESeq_expression_data_controlIsUnreated.csv', x=joint)
write.csv(file='DESeq_norm_size_factor_data.csv', x=normalized)
cat('Saving results in DESeq_expression_data_controlIsUnreated.csv. ', '\n')

# Identification of significantly expressed mi/pi/sno-RNA
cat('DESeq2 analysis is complete. Now finding significantly expressed entries.....','\n')
input_mir <- read.csv('DESeq_expression_data_controlIsUnreated.csv')

# Selecting differentially expressed miRNAs with p_adj <= 0.05
output_mir <- (subset(input_mir,input_mir$padj <= 0.05,select=cts...1.))

# Saving results
write.table(file='significantly_DE_mir.csv', x=output_mir,quote=FALSE,row.names=FALSE,col.names=FALSE)
cat('Final the results are saved in significantly_DE_mir.csv ')

