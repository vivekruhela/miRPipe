cat('Checking if BiocManager is install...')
if (!requireNamespace("BiocManager", quietly = TRUE)) { 
	install.packages("BiocManager")
    BiocManager::install(c("DESeq2","ShortRead","Biostrings"))
} else {
	cat('BiocManager and other packages are installed')
}
