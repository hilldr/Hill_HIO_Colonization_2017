## Differential expression of kallisto results with DESeq2

## Retrieve kallisto counts for ECOR2 0-96 hr vs PBS dataset (Figure 3)
## read in table with sample metadata
samples <- readr::read_csv(file = "../data/RNA-seq/sample_key.csv")

## subsetting rules
samples <- subset(samples, (samples$Date == "2016-08-18" &
                                ((samples$injection == "Escherichia coli str. ECOR2" & 
                                                                              samples$heat_killed == 0 &
                                                                                  samples$PMN == 0 &
                                                                                  samples$hr == 24 &
                                                                                      samples$hypoxia == 0 &
                                                                                      samples$pharmacologic == "none") |
                                 (samples$injection == "PBS" & 
                                      samples$heat_killed == 0 &
                                      samples$PMN == 0 &
                                          samples$hr == 24 &
                                          samples$hypoxia == 0 &
                                              samples$pharmacologic == "none"))) |
                               (samples$hr == 48 & samples$injection != "PBS") | samples$hr == 96)

## set PBS to 0 hr                  
samples[samples$injection == "PBS",]$hr <- 0

## setup access to kallisto read files
files <- file.path(samples$directory,
                   samples$file_name,
                   "abundance.tsv")

## set sample names as description_rep#_seq_rep#
names(files) <- paste0(samples$code_name,"-",samples$hr,"_",samples$num)
## check that all files are found
if (all(file.exists(files)) == FALSE) {
    print("kallisto files not found")
    stop()
}

## associate transcripts with gene IDs
## check if saved transcript:gene index is present
## recommended - biomaRt connectivity is unreliable
if (file.exists("../data/RNA-seq/tx2gene.Rdata") == TRUE) {
    load(file = "../data/RNA-seq/tx2gene.Rdata")
} else {    
    ## create biomart reference
    ensembl <- biomaRt::useMart("ensembl")
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'useast.ensembl.org')
    ## create index of gene names
    tx2gene <- biomaRt::getBM(attributes = c("ensembl_transcript_id","external_gene_name"), mart = mart)
}


## import kallisto data and generate count dataframe (dds)
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
library(readr)
txi <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)

## create directory to deposit results
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"
dir.create(path = data.dir, recursive = TRUE)
## export transcript abundance counts
write.csv(txi$abundance, file = file.path(data.dir,"complete_dataset_txi.csv"))

library(DESeq2)
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
dds <- DESeq2::DESeqDataSetFromTximport(txi,
                                        colData = samples,
                                        design = ~ hr) #code_name
## pre-filter out counts < 1
dds <- dds[rowSums(counts(dds)) > 0.1, ]

## write out normalized expression counts
dds <- DESeq2::estimateSizeFactors(dds)
ddscounts <- DESeq2::counts(dds, normalized = TRUE)

## write expression matrix to file
write.csv(ddscounts,
          file =  file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv"))

## Generate differential expression results in DESeq2

## enable parallel processes
#library("BiocParallel")
#register(MulticoreParam(4))

data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"

colData(dds)@listData$hr[colData(dds)@listData$injection == "PBS"] <- 0

## setup multifactor design
colData(dds)$code_name <- as.factor(paste0(colData(dds)$code_name, "_", colData(dds)$hr))
ddsMF <- dds
DESeq2::design(ddsMF) <- ~ code_name

## Likelihood ratio test (ANOVA-like)
## set parallel
ddsMF <- DESeq2::DESeq(ddsMF, test = "LRT", reduced = ~1, parallel = FALSE)
res <- DESeq2::results(ddsMF)
write.csv(res, file = file.path(data.dir, "LRT.csv"))

## Wald tests
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "ECOR2_24", "PBS_0"))
write.csv(res, file = file.path(data.dir, "ECOR2_over_PBS_24hr.csv"))

res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "ECOR2_48", "PBS_0"))
write.csv(res, file = file.path(data.dir, "ECOR2_over_PBS_48hr.csv"))

res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "ECOR2_96", "PBS_0"))
write.csv(res, file = file.path(data.dir, "ECOR2_over_PBS_96hr.csv"))
