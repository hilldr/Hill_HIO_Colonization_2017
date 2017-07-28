## Figure 6 multipanel ---------------------------------------------------------
source("ggplot2-themes.R")
source("custom_fun.R")
library(ggplot2)
library(gridExtra)

sfig4a <- png2ggplot("../data/figure7/occludin_cropped.png") +
    img.theme + coord_fixed(ratio = 1/1.6) + ggtitle("A")

sfig4b <- png2ggplot("../data/figure7/acTub_cropped.png") +
    img.theme + coord_fixed(ratio = 1/1.6) + ggtitle("B")


layout <- rbind(c(1),
                c(2))

## PDF output
pdf(file = "../figures/supplemental/sfigure4_multipanel.pdf", width = 4125/300, height = 5250/300, onefile = FALSE)
gridExtra::grid.arrange(sfig4a,sfig4b, layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/supplemental/sfigure4_multipanel.eps", 
       plot = gridExtra::grid.arrange(sfig4a,sfig4b, layout_matrix = layout), 
       width = 22, height = 28)

## PNG output
png(filename = "../figures/supplemental/sfigure4_multipanel.png", width = 1100, height = 1400)
gridExtra::grid.arrange(sfig4a,sfig4b, layout_matrix = layout)
dev.off()

## Supplemental FIGURE  --------------------------------------------------------------------
## supplemental figure 1A
## import data
## Differential expression of kallisto results with DESeq2

## read in table with sample metadata
samples <- readr::read_csv(file = "../data/RNA-seq/sfigure1_sample_sheet.csv")

## setup access to kallisto read files
files <- file.path(samples$directory,
                   samples$file_name,
                   "abundance.h5")

## set sample names as description_rep#_seq_rep#
names(files) <- paste0(samples$description,"_",
                       samples$rep,"_L00",samples$seq_rep)

## check that all files are found
if (all(file.exists(files)) == FALSE) {
    print("kallisto files not found")
    stop()
}

## associate transcripts with gene IDs
## create biomart reference
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')
## create index of gene names
tx2gene <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                         "external_gene_name"),
                          mart = mart)

## import kallisto data and generate count dataframe (dds)
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
library(readr)
txi <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene)

data.dir <- "../results/supplemental"
## export abundance counts
write.csv(txi$abundance, file = file.path(data.dir, "tximport_abundance counts.csv"))

library(DESeq2)
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
dds <- DESeq2::DESeqDataSetFromTximport(txi,
                                        colData = samples,
                                        design = ~ status)
## collapse technical replicates
ddsColl <- DESeq2::collapseReplicates(object = dds,
                                      groupby = dds$sample_id, 
                                      run = dds$seq_rep,
                                      renameCols = TRUE)

## pre-filter out counts < 1
dds <- dds[rowSums(counts(dds)) > 1, ]
ddsColl <- ddsColl[rowSums(counts(ddsColl)) > 1, ]

## write out normalized expression counts
ddsColl <- DESeq2::estimateSizeFactors(ddsColl)
ddsCollcounts <- DESeq2::counts(ddsColl, normalized = TRUE)
## give the columns better names
colnames(ddsCollcounts) <- paste0(colData(ddsColl)@listData$status,"_",
                                  colData(ddsColl)@listData$origin,"_",
                                  colData(ddsColl)@listData$rep)
## write expression matrix to file
write.csv(ddsCollcounts, file =  file.path(data.dir, "DESeq2-normalized-counts.csv"))
save(dds, file = file.path(data.dir, "dds.Rdata"))
rm(list = ls()) 

## Pearson's correlation matrix ------------------------------------------------
## Load dataset from file
data.dir <- "../results/supplemental"

df <- readr::read_csv(file.path(data.dir, "DESeq2-normalized-counts.csv"))


num.data <- df[,sapply(df,is.numeric)]
group <- gsub('.{2}$', '', colnames(num.data))

## pca.data <- num.data[apply(num.data, 1, sd, na.rm=TRUE) != 0,]
## calculate variance by row (gene)
var <- apply(num.data, 1, sd, na.rm=TRUE)
## adjust cut off according to variance percentile
pca.data <- num.data[var > quantile(var, 0.2) & var != 0,]
colnames(pca.data) <- gsub("hES_hES", "hPSC", colnames(pca.data))
colnames(pca.data) <- gsub("_", " ", colnames(pca.data))
colnames(pca.data) <- gsub("in vivo", "tx", colnames(pca.data))
cor1 <- cor(pca.data, method = "pearson")

## determine order for axis clustering
library(magrittr)
library(ggtree)
library(ape)
tree <- dist(cor1, method = "canberra") %>% hclust(method = "mcquitty") %>% as.phylo
plot <- ggplot(tree) + geom_tree(size = 3) + theme_tree() + geom_tiplab(size = 4)

library(gtable)
library(ggplot2)
library(grid)
source("ggplot2-themes.R")
## Correlation matrix
ord <- hclust(dist(cor1, method = "canberra"), method = "mcquitty")$order
melted_cormat <- reshape2::melt(cor1[ord,ord])
sfigure1a <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
##  geom_tile(color = "grey", size = 0.5) +  # add tiles?
    geom_tile() +
    scale_fill_distiller(expression(paste(italic(r))),palette = "Spectral") +
    xlab("") + ylab("") + coord_fixed(ratio = 1) + theme1 +
    theme(axis.text = element_text(size = 18, face ="bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.title = element_text(size = 18, face ="bold"),
	  panel.background = element_rect(fill = "white"),
	  panel.border = element_blank(),
          legend.text = element_text(size = 12, face ="bold")) + ggtitle("A") 
## tree <- dist(cor1, method = "canberra") %>% hclust(method = "mcquitty") %>% as.phylo
## tree.plot <- ggplot(tree) + geom_tree(size = 2) + theme_tree() + geom_tiplab(size = 0)

## ## left strip of labels
## b.labels <- data.frame(names = rownames(cor1)[ord], x = c(1:18))
## p <- ggplot(b.labels, aes(xmin=0, xmax=1, ymin=1, ymax=x+1)) + geom_rect(fill = "white")
## p <- p + geom_text(aes(y=x+0.5, x=1, label=paste0(names, " ")), hjust=1, vjust=0.5, angle=0, size = 5, fontface = "bold") 
## p <- p + scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
## pleft <- p

## ## setp grid plot
## g <- gtable_filter(ggplotGrob(plot), pattern = "panel", trim = TRUE, fixed=TRUE)
## g <- gtable_add_cols(g, unit(0.25, "null"), 0) # left labels
## g <- gtable_add_cols(g, unit(0.2, "null"), 1) # left labels
## g <- gtable_add_cols(g, unit(0.15, "null"), -1) # left labels
## g <- gtable_add_rows(g, unit(0.15, "null"), 4) # bottom labels
## g <- gtable_add_grob(g, gtable_filter(ggplotGrob(tree.plot), pattern = "panel", trim = TRUE, fixed=TRUE), 1,1)
## g <- gtable_add_grob(g, gtable_filter(ggplotGrob(pleft), pattern = "panel", trim = TRUE, fixed=TRUE), 1,2)
## g <- gtable_add_grob(g, gtable_filter(ggplotGrob(plot), pattern = "axis-b", trim = TRUE, fixed=TRUE), 2,3)
## g <- gtable_add_grob(g, gtable_filter(ggplotGrob(plot), pattern = "guide-box", trim = TRUE, fixed=TRUE), 1,4)
#grid.newpage()
#grid.draw(g)
 
png(filename = "../figures/supplemental/sfigure1a.png", width = 900, height = 900)
print(sfigure1a)
 dev.off()
 ## ggsave(filename = "../figures/supplmental/eps/sfigure1a.eps", 
 ## plot = grid.draw(g), 
 ## width = 20, height = 20)

library(magrittr)
cordf <- cor1 %>% as.data.frame() %>% dplyr::select(dplyr::matches("HIO|fetal|adult"))
## select adult SI comparisons
set <- grep("adult", rownames(cordf))
cordf <- cordf[set,]
## reshape for easy ggplot
cordf <- cordf %>% reshape2::melt()
## generic group names
cordf$variable <- gsub('.{2}$', '', cordf$variable)
cordf$variable <- factor(cordf$variable,
                         levels = c("in vitro HIO", "fetal SI","tx HIO", "adult SI", "hPSC"),
                         ordered = TRUE) 

## stats
test <- t.test(x = cordf[cordf$variable == "in vitro HIO",]$value,
               y= cordf[cordf$variable == "tx HIO",]$value,
               alternative = "two.sided")

cordf2 <- cor1 %>% as.data.frame() %>% dplyr::select(dplyr::matches("HIO|fetal|adult|hPSC"))
set <- grep("fetal SI", rownames(cordf2))
cordf2 <- cordf2[set,]
## reshape for easy ggplot
cordf2 <- cordf2 %>% reshape2::melt()
## generic group names
cordf2$variable <- gsub('.{2}$', '', cordf2$variable)
cordf2$variable <- factor(cordf2$variable,
                         levels = c("in vitro HIO", "fetal SI","tx HIO", "adult SI", "hPSC"),
                         ordered = TRUE)

## These tests are referenced in the text
# format(test2$p.value, digits = 3)
test2 <- t.test(x = cordf2[cordf2$variable == "in vitro HIO",]$value,
               y= cordf2[cordf2$variable == "hPSC",]$value,
               alternative = "two.sided") 

test3 <- t.test(x = cordf2[cordf2$variable == "tx HIO",]$value,
               y= cordf2[cordf2$variable == "hPSC",]$value,
               alternative = "two.sided")

library(ggplot2)
source("ggplot2-themes.R")

sfigure1b <- ggplot(data = cordf, aes(x = variable, y = value, color = variable)) +
    geom_boxplot(size = 2, outlier.size = 5) + theme1 +
    ylab(expression(paste("Pearson's correlation coefficent (", italic(r),")"))) +
    xlab("Comparison vs. adult SI") +
    scale_color_manual(values = c(color.set[2], color.set[3], color.set[4], color.set[1])) + theme1 +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    annotate("segment", x = 1, xend = 3, y = 0.96, yend = 0.96, color = "black", size = 2) +
    annotate("text", x = 2, y = 0.965, size = 10, label = paste("P =", round(test$p.value, 5))) +
    ggtitle("B")

png(filename = "../figures/supplemental/sfigure1b.png", width = 500, height = 1000)
print(sfigure1b)
dev.off()

## PCA analysis ------------------------------------------------------
## Load dataset from file
data.dir <- "../results/supplemental"
df <- readr::read_csv(file.path(data.dir, "DESeq2-normalized-counts.csv"))
## rename for SeqRetriever compatibility
df <- df[,-grep("hES", colnames(df))]
num.data <- df[,sapply(df,is.numeric)]
group <- gsub('.{2}$', '', colnames(num.data))
group <- gsub("_", " ", group)
group <- gsub("in vivo", "tx", group)

## pca.data <- num.data[apply(num.data, 1, sd, na.rm=TRUE) != 0,]
## calculate variance by row (gene)
var <- apply(num.data, 1, sd, na.rm=TRUE)
## adjust cut off according to variance percentile
pca.data <- num.data[var > quantile(var, 0.2) & var != 0,]
cor1 <- cor(pca.data, method = "pearson")
pca <- prcomp(t(pca.data),scale = TRUE,center = TRUE)
scores <- data.frame(colnames(pca.data), pca$x[,1:ncol(pca$x)],group)

index <- data.frame(group = as.vector(unique(scores$group)),
                    ellipse = c("mature", "immature","immature", "mature"))

scores <- plyr::join(scores, index, by = 'group')


## PCA plot
## function to format decimals as precentage
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

library(ggplot2)
library(RColorBrewer)
sfigure1c <- qplot(x = PC1, y = PC2, data = scores) +  
    geom_point(shape = 21,aes(fill = factor(group)), size = 10, color = "white") +
    scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(scores$group)), "Set1"))(length(unique(scores$group)))) +
    theme1 +
   # theme(legend.position="bottom",
         # legend.background = element_rect(colour = "white"),
        #  legend.key = element_rect(color = "white",fill = "white"),
#	  panel.border = element_rect(fill = NA, color = "white")) +
    theme(legend.position="bottom",
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(color = "white",fill = "white"),
	  panel.border = element_rect(fill = NA, color = "grey70"),
          legend.title = element_text(size = 32),
          legend.text = element_text(size = 24)) +
    guides(fill = guide_legend(title = "Epithelium", nrow = 4, byrow=TRUE, order =1),
           color = guide_legend(title = "Status", nrow = 2, byrow=TRUE, order =2)) +
    coord_fixed(ratio = 1) +
    xlab(paste("PC1 (",percent(round(summary(pca)$importance[2,1],4)),")",sep = "")) +
    ylab(paste("PC2 (",percent(round(summary(pca)$importance[2,2],4)),")",sep = "")) +
    stat_ellipse(type = "t", linetype = 2, size = 2, aes(color = factor(ellipse))) + 
    ggtitle("C")

png(filename = "../figures/supplemental/sfigure1c.png", width = 1000, height = 1000)
print(sfigure1c)
dev.off()

## Enteroids ONLY processing
## Differential expression of kallisto results with DESeq2

## read in table with sample metadata
## read in table with sample metadata
samples <- readr::read_csv(file = "../data/RNA-seq/sfigure1_sample_sheet.csv")
## remove hES samples
samples <- subset(samples, samples$description != "hES")

## setup access to kallisto read files
files <- file.path(samples$directory,
                   samples$file_name,
                   "abundance.h5")

## set sample names as description_rep#_seq_rep#
names(files) <- paste0(samples$description,"_",
                       samples$rep,"_L00",samples$seq_rep)

## check that all files are found
if (all(file.exists(files)) == FALSE) {
    print("kallisto files not found")
    stop()
}

## associate transcripts with gene IDs
## create biomart reference
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')
## create index of gene names
tx2gene <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                         "external_gene_name"),
                          mart = mart)

## import kallisto data and generate count dataframe (dds)
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
library(readr)
txi <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene)

## export abundance counts
data.dir <- "../results/supplemental"
write.csv(txi$abundance, file = file.path(data.dir,"tximport_abundance counts_noPSC.csv"))

library(DESeq2)
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
dds <- DESeq2::DESeqDataSetFromTximport(txi,
                                        colData = samples,
                                        design = ~ status)
## collapse technical replicates
ddsColl <- DESeq2::collapseReplicates(object = dds,
                                      groupby = dds$sample_id, 
                                      run = dds$seq_rep,
                                      renameCols = TRUE)

## pre-filter out counts < 1
dds <- dds[rowSums(counts(dds)) > 1, ]
ddsColl <- ddsColl[rowSums(counts(ddsColl)) > 1, ]

## write out normalized expression counts
ddsColl <- DESeq2::estimateSizeFactors(ddsColl)
ddsCollcounts <- DESeq2::counts(ddsColl, normalized = TRUE)
## give the columns better names
colnames(ddsCollcounts) <- paste0(colData(ddsColl)@listData$status,"_",
                                  colData(ddsColl)@listData$origin,"_",
                                  colData(ddsColl)@listData$rep)
## write expression matrix to file
write.csv(ddsCollcounts, file =  file.path(data.dir,"DESeq2-normalized-counts_noPSC.csv"))
save(ddsColl, file = file.path(data.dir, "dds2.Rdata"))
rm(list = ls())

## enable parallel processes
library("BiocParallel")
register(MulticoreParam(4))
load("../results/supplemental/dds2.Rdata")
## setup multifactor design
ddsMF <- ddsColl
DESeq2::design(ddsMF) <- ~ status

## Likelihood ratio test (ANOVA-like)
ddsMF <- DESeq2::DESeq(ddsMF, test = "LRT", reduced = ~1, parallel = TRUE)
res <- DESeq2::results(ddsMF)

## need to specify Wald test when later making specific comparisons
res.ivit <- DESeq2::results(ddsMF, test ="Wald",
                           contrast = c("status", "in_vivo", "in_vitro"))
res.af <- DESeq2::results(ddsMF, test ="Wald",
                           contrast = c("status", "adult", "fetal"))
## write out stats results
data.dir <- "../results/supplemental"
write.csv(res.ivit, file = file.path(data.dir,"in-vivo_in-vitro_Wald.csv"))
write.csv(res.af, file = file.path(data.dir,"adult_fetal_Wald.csv"))

## Expression plot -------------------------------------------------------------
## Load dataset from file
data.dir <- "../results/supplemental"
## adult/fetal
adfe <- readr::read_csv(file = file.path(data.dir,"in-vivo_in-vitro_Wald.csv"))
## in vivo/in vitro
ivit <- readr::read_csv(file = file.path(data.dir,"adult_fetal_Wald.csv"))

## subset columns of interest
library(magrittr)
adfe <- adfe %>% dplyr::select(X1, log2FoldChange, padj) %>% dplyr::rename(adfe_log2 = log2FoldChange, adfe_padj = padj)
ivit <- ivit %>% dplyr::select(X1, log2FoldChange, padj) %>% dplyr::rename(ivit_log2 = log2FoldChange, ivit_padj = padj)

## join dataframes
df <- dplyr::left_join(adfe, ivit, by = 'X1')

## subset data prior to plotting
df <- subset(df, df$adfe_padj < 0.05 | df$ivit_padj < 0.05, na.rm = TRUE)

## number of DE genes on each axis, reference in text
num.DE.si <- length(rownames(subset(df, df$adfe_padj < 0.05)))
num.DE.hio <- length(rownames(subset(df, df$ivit_padj < 0.05)))
## number of fetal signature genes
num.fetal <- length(rownames(subset(df, df$adfe_log2 < 0 &  df$ivit_log2 < 0)))
## number of adult signature genes
num.adult <- length(rownames(subset(df, df$adfe_log2 > 0 &  df$ivit_log2 > 0)))

## create adult/in vitro gene set
adult.tHIO.gs <- df[df$adfe_log2 > 0 &
                        df$ivit_log2 > 0 &
                        df$adfe_padj < 0.001  &
                        df$ivit_padj < 0.001,] %>% dplyr::rename(SYMBOL = X1)
write.csv(adult.tHIO.gs, file = file.path(data.dir,"adult-tHIO-gene-set.csv"))

## create adult/in vitro gene set
adult.tHIO.gs <- df[df$adfe_log2 > 1 &
                        df$ivit_log2 > 1 &
                        df$adfe_padj < 0.05  &
                        df$ivit_padj < 0.05,] %>% dplyr::rename(SYMBOL = X1)
write.csv(adult.tHIO.gs, file = file.path(data.dir,"adult-tHIO-gene-set2.csv"))

## create adult/in vitro gene set
adult.tHIO.gs <- df[df$adfe_log2 > 0 &
                        df$ivit_log2 > 0 &
                        df$adfe_padj < 0.05  &
                        df$ivit_padj < 0.05,] %>% dplyr::rename(SYMBOL = X1)
write.csv(adult.tHIO.gs, file = file.path(data.dir,"adult-tHIO-gene-set3.csv"))

## linear regression
model1d <- lm(df$adfe_log2 ~ df$ivit_log2)

library(ggplot2)
source("ggplot2-themes.R")

sfigure1d <- ggplot(data = df, aes(x = ivit_log2, y = adfe_log2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", na.rm = TRUE, contour = FALSE) +
    xlim(c(-1.5, 1.5)) +
    ylim(c(-1.5,1.5)) +
    theme1 +
    scale_fill_distiller("Density", palette = "Blues", direction = 1) +
    labs(x = expression(Log[2]*"FC(in vivo HIO / in vitro HIO)"),
         y = expression(Log[2]*"FC(adult SI / fetal SI)")) +
    annotate("segment", size =1, x = 0, xend = 0,  y = -1.5, yend =1.5, color = "white", linetype = "dashed") +
    annotate("segment", size =1, x = -1.5, xend = 1.5,  y = 0, yend = 0, color = "white", linetype = "dashed") +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 32),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 12),
          panel.border = element_rect(fill = NA, color = "white")) +
    ggtitle("D")

png(filename = "../figures/supplemental/sfigure1d.png", width = 1000, height = 1000)
print(sfigure1d)
dev.off()

## Figure 1 multipanel ---------------------------------------------------------
library(ggplot2)

## PDF output
pdf(file = "../figures/supplemental/sfigure1_multipanel.pdf", width = 8400/300, height = 8400/300, onefile = FALSE)
gridExtra::grid.arrange(sfigure1a, sfigure1b, sfigure1c, sfigure1d, 
             layout_matrix = rbind(c(1,1,2,2),
                                   c(3,3,4,4)))
dev.off()

## PNG output
png(filename = "../figures/supplemental/sfigure1_multipanel.png", width = 2100, height = 2100)
gridExtra::grid.arrange(sfigure1a, sfigure1b, sfigure1c, sfigure1d, 
             layout_matrix = rbind(c(1,1,2,2),
                                   c(3,3,4,4)))
dev.off()
