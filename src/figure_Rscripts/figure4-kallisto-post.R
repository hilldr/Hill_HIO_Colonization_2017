## Differential expression of kallisto results with DESeq2

## read in table with sample metadata
samples <- readr::read_csv(file = "../data/RNA-seq/sample_key.csv")

samples <- subset(samples, samples$Date == "2016-08-18" | (samples$Date == "2016-05-24" & samples$injection != "hypoxia"))

write.csv(samples, "../data/RNA-seq/figure4_sample_key.csv")

## setup access to kallisto read files
files <- file.path(samples$directory,
                   samples$file_name,
                   "abundance.tsv")

## set sample names as description_rep#_seq_rep#
names(files) <- samples$short_name
## check that all files are found
if (all(file.exists(files)) == FALSE) {
    print("kallisto files not found")
    stop()
}

## associate transcripts with gene IDs
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
data.dir <- "../results/ECOR2_hypoxia_nfkb/"
dir.create(path = data.dir, recursive = TRUE)
## export abundance counts
write.csv(txi$abundance, file = file.path(data.dir,"complete_dataset_txi.csv"))

library(DESeq2)
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
dds <- DESeq2::DESeqDataSetFromTximport(txi,
                                        colData = samples,
                                        design = ~ code_name)
## pre-filter out counts < 1
dds <- dds[rowSums(counts(dds)) > 0.1, ]

## write out normalized expression counts
dds <- DESeq2::estimateSizeFactors(dds)
ddscounts <- DESeq2::counts(dds, normalized = TRUE)

## write expression matrix to file
write.csv(ddscounts,
          file =  file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv"))

## enable parallel processes
library("BiocParallel")
register(MulticoreParam(4))

## setup multifactor design
colData(dds)$code_name <- as.factor(colData(dds)$code_name)
ddsMF <- dds
DESeq2::design(ddsMF) <- ~ code_name

## Likelihood ratio test (ANOVA-like)
## check parallel = FALSE if error occurs at this point (this was broken by upgrade to R 3.3.3)
ddsMF <- DESeq2::DESeq(ddsMF, test = "LRT", reduced = ~1, parallel = TRUE)
res <- DESeq2::results(ddsMF)
write.csv(res, file = file.path(data.dir, "LRT.csv"))

## need to specify Wald test when later making specific comparisons
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "ECOR2", "ECOR2+NFKB-inhib"))
write.csv(res, file = file.path(data.dir, "ECOR2_over_ECOR-NFKBi.csv"))
## vs pbs
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "ECOR2", "PBS"))
write.csv(res, file = file.path(data.dir, "ECOR2_over_PBS.csv"))
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "ECOR2-heat-kill", "PBS"))
write.csv(res, file = file.path(data.dir, "ECOR2-HK_over_PBS.csv"))
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "hypoxia", "PBS"))
write.csv(res, file = file.path(data.dir, "hypoxia_over_PBS.csv"))
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "ECOR2+NFKB-inhib", "PBS"))
write.csv(res, file = file.path(data.dir, "ECOR2-NFKBi_over_PBS.csv"))
## vs. Ecor2
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "PBS", "ECOR2"))
write.csv(res, file = file.path(data.dir, "PBS_over_ECOR2.csv"))
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "ECOR2-heat-kill", "ECOR2"))
write.csv(res, file = file.path(data.dir, "ECOR2-HK_over_ECOR2.csv"))
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "hypoxia", "ECOR2"))
write.csv(res, file = file.path(data.dir, "hypoxia_over_ECOR2.csv"))

res <- DESeq2::results(ddsMF, test = "Wald",
                        contrast = c("code_name", "hypoxia", "hypoxia+NFKB-inhib"))
write.csv(res, file = file.path(data.dir, "hypoxia_over_hypoxia-NFKBi.csv"))

res <- DESeq2::results(ddsMF, test = "Wald",
                        contrast = c("code_name", "ECOR2-heat-kill", "ECOR2-heat-kill+NFKB-inhib"))
write.csv(res, file = file.path(data.dir, "ECOR2-HK_over_ECOR2-HK-NFKBi.csv"))
