## Figure 8 multipanel ---------------------------------------------------------
source("ggplot2-themes.R")
source("custom_fun.R")
library(ggplot2)
library(gridExtra)

sfig4a <- png2ggplot("../data/figure8/occludin_cropped.png") +
    img.theme + coord_fixed(ratio = 1/1.6) + ggtitle("A")

sfig4b <- png2ggplot("../data/figure8/acTub_cropped.png") +
    img.theme + coord_fixed(ratio = 1/1.6) + ggtitle("B")


layout <- rbind(c(1),
                c(2))

## PDF output
pdf(file = "../figures/supplemental/sfigure5_multipanel.pdf", width = 4125/300, height = 5250/300, onefile = FALSE)
gridExtra::grid.arrange(sfig4a,sfig4b, layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/supplemental/sfigure5_multipanel.eps", 
       plot = gridExtra::grid.arrange(sfig4a,sfig4b, layout_matrix = layout), 
       width = 22, height = 28)

## PNG output
png(filename = "../figures/supplemental/sfigure5_multipanel.png", width = 1100, height = 1400)
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

library(ggplot2)
library(ggtree)
library(magrittr)
source("ggplot2-themes.R")
tree <- read.tree(file = "../data/ecor2/ECOR2-50.nwk")
plot <- ggtree(tree, size = 2) %>% ggtree::collapse(node = 86)
plot <- plot + geom_treescale(width = 0.015) + #xlim(NA, 100) +
  #  geom_hilight(76, "steelblue") +
    geom_tiplab(size = 7) +
    geom_cladelabel(node = 53,
                    label = "non-pathogenic E. coli",
                    align = FALSE,
                    offset = 0.01,
                    offset.text = 1e-3,
                    fontsize = 12,
                    barsize = 2) +
    #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
    geom_tippoint(aes(subset=(node == 26)),
                  size = 8,
                  shape = 21,
                  fill = color.set[1],
                  color = color.set[1])

    
## PDF output
pdf(file = "../figures/figure1/sfigure1-2_tree.pdf", width = 9000/300, height = 7500/300, onefile = FALSE)
print(plot)
dev.off()

## Figure 1 - Supplement 3
# Import dataset in SRC header
library(ggplot2)
source("ggplot2-themes.R")
library(scales)      

cfudata <- readr::read_csv(file = "../data/figure1/culture_optimization/LB_pbs_data.csv")

# T-test
stats <- t.test(cfudata[cfudata$media == "LB",]$CFU_HIO,cfudata[cfudata$media == "PBS",]$CFU_HIO,alternative= "greater", var.equal =TRUE)
lbplot <- ggplot(cfudata, aes( x= factor(media), y = CFU_HIO)) +
    geom_boxplot(aes(color = factor(media)), size = 2) +
   # geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    ylab("CFU/HIO") + xlab("") +
    annotate("text", y=4000000, x= 2,
             label = paste("P = ", format(stats$p.value,digits = 1)),
             color="grey30",size =10) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme1 + ggtitle("B")

print(lbplot)

data <- read.table("../data/figure1/culture_optimization/ECOR2growth.csv",
                   header = TRUE, sep = ",",
                   stringsAsFactors=FALSE)
stats <- t.test(data[data$sample == "Remove abx",]$mean,
                data[data$sample == "Retain abx",]$mean,
                alternative= "greater", var.equal =TRUE)
library(ggplot2)
source("ggplot2-themes.R")
library(scales)      
fig1b <- ggplot(data, aes(x=sample, y=mean,)) +
    geom_boxplot(aes(color=sample), size = 2) +
    ylab("CFU/HIO") +
    xlab("") + theme1 + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="l",size= 2) +
    annotate("text", y=400000, x= 2,
             label = paste(" P = ",round(stats$p.value,digits = 3)),
             color = "grey30", size = 8) + ggtitle("C")

print(fig1b)

## import data
data <- readr::read_csv(file = "../data/figure1/culture_optimization/zoi_data.csv", col_names = TRUE)

## convert measured pixels to um
## 1mm = 106 px
data$mm <- data$Length/106

## add group designations
data$treatment <- rep(c("ENR + abx", "ENR after washout"), times = length(data$Label)/2)

library(ggplot2)
source("ggplot2-themes.R")
source("custom_fun.R")

library(ggstance)

data$treatment <- factor(data$treatment, levels = c( "ENR after washout", "ENR + abx"))

stats <- t.test(data[data$treatment == "ENR after washout",]$mm,
                data[data$treatment == "ENR + abx",]$mm,
                alternative= "less", var.equal =TRUE)

zoi <- ggplot(data, aes(y = mm, x = treatment,)) +
    geom_boxplot(aes(color=treatment), size = 1) +
    xlab("") +
    ylab("Inhibition zone diameter (mm)") + theme1 + 
    annotation_logticks(sides="l",size= 2) +
    annotate("text", y= 1, x= 1,
             label = paste("P = ", format(stats$p.value,digits = 1)),
             color="grey30",size =10) +
    theme1 +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust =1))

fig <- png2ggplot("../data/figure1/culture_optimization/representative_zoi.png") +
    img.theme + coord_fixed(ratio = 1/2) + ggtitle("D")

layout2 <- rbind(c(1),
                c(2),
                c(2))


gridExtra::grid.arrange(fig, zoi, layout_matrix = layout2)

## Figure 1 - Supplement 3 multipanel ---------------------------------------------------------
library(ggplot2)
library(gridExtra)
source("ggplot2-themes.R")
source("custom_fun.R")

schematic <- png2ggplot(filename = "../figures/figure1/schematic.png") +
    coord_fixed(ratio = 1.5) + ggtitle("A") + img.theme

layout <- rbind(c(1,1,1,2),
                c(1,1,1,3),
                c(1,1,1,4),
                c(1,1,1,4))


sfig3d <- gridExtra::grid.arrange(fig, zoi, layout_matrix = layout2)

## PDF output
pdf(file = "../figures/figure1/sfigure1-3_multipanel.pdf", width = 7500/300, height = 7500/300, onefile = FALSE)
gridExtra::grid.arrange(schematic, lbplot, fig1b, sfig3d, 
             layout_matrix = layout)
dev.off()

#+begin_src R :session *R* :exports none :results graphics :file ../figures/supplemental/sfigure1_supp2.pdf :eval yes :tangle figure_Rscripts/supplemental.R
## Minimum density /E. coli/ required to establish colonization 

## E. coli growth in HIOs
## Read in data table
data <- read.table("../data/figure1/ECOR2growth_fig1.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
## Index of CFU/HIO injected for each Sample condition (A-G)
sample.table <- read.table("../data/figure1/sample_table_fig1.csv", header = TRUE, sep = ",", stringsAsFactors=FALSE)
## Generate index of rows in sample table that match
## the sample labels in data
id <- match(data$sample,sample.table$sample)
## create column in data of of CFU/HIO values in sample table
## in matching rows listed in id
data$inject <- sample.table[id,]$value
data$fold <- data$mean/data$inject
data$increase <- ifelse(data$mean > data$inject,"increase","decrease")
  
group <- aggregate(mean ~ inject, data = data, FUN = mean)
group.sem <- aggregate(mean ~ inject, data = data,
                       FUN = function(x) sd(x)/sqrt(length(x)))
group$sem <- group.sem$mean

## ANOVA of mean CFU/HIO among colonized HIOs
fit <- aov(mean ~ sample, data = data[data$inject >1,])
fit2 <- lm(log(data[data$inject > 0 & data$fold !=0,]$fold) ~ log(log(data[data$inject > 0 & data$fold !=0,]$inject)), data[data$inject > 0 & data$fold !=0,])

p1c <- summary(fit)[[1]][["Pr(>F)"]][[1]]

## % colonized at < 5 CFU
pct1 <- round(100*(1-length(rownames(data[data$inject < 5 & data$inject > 1 & data$mean < 1,] ))/length(rownames(data[data$inject < 5 & data$inject > 1 & data$mean > 1,] ))),1)

## % colonized at > 100 CFU
pct2 <- round(100*(1-length(rownames(data[data$inject > 100 & data$mean < 1,] ))/length(rownames(data[data$inject < 5 & data$inject > 1 & data$mean > 1,] ))),2)

## plot
library(ggplot2)
library(grid)
library(scales)
source("ggplot2-themes.R")

scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

## generate stats string for plot
source("custom_fun.R")
stats <- lm_eqn(data[data$inject > 0 & data$fold !=0,],
                data[data$inject > 0 & data$fold !=0,]$inject,
                data[data$inject > 0 & data$fold !=0,]$mean)


fig2c <- ggplot(data, aes(x=inject, y=mean)) +
    geom_smooth(
                aes(x=inject, y=mean), colour = "black",
                size = 2,
                method = "lm",
                formula = y ~ x,
                level = 0.95) +
  #  geom_boxplot(aes(group = inject), size =2, fill = color.set[2]) +
    geom_point(size = 8, shape =21, fill = color.set[2]) +
    scale_fill_brewer(palette = "Set1") + 
    ylab(latex2exp::TeX("$\\textbf{CFU$\\cdot{}HIO^{-1}_{$\\textit{t}=24}}$")) +
    xlab("CFU injected per HIO") + theme1 +
    scale_y_log10(limits = c(1,50000000),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(limits = c(1,100000),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +    
    annotation_logticks(sides = "bl", size = 2,
                        short = unit(.75,"mm"),
                        mid = unit(3,"mm"),
                        long = unit(5,"mm")) +
        ## size of stats label
    annotate("text",x = 5e3, y = 5e7,
                           label = substr(stats,62,150),
                           parse = TRUE,
                           size = 10) +


pdf(file = "../figures/supplemental/sfigure1_supp2.pdf", width = 10, height = 10, onefile = FALSE)
print(fig2c)
dev.off()

## FIGURE 5 - Supplement 1 ---------------------------------------------------
## Figure 5 Supplement 1A: RNA-seq NF-kb heatmap
## import data
## Load dataset from file
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"
library(magrittr)
df <- readr::read_csv(file = file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv")) %>% dplyr::rename(SYMBOL = X1)

## list of genes to extract
paths <- readr::read_csv("../results/ECOR2_hypoxia_nfkb/plotted-nfkb_complete-goANDreactome-results.csv")

nfkb <- paths[grep("NF-kB", paths$Description),]$geneID
nfkb.path <- paths[grep("NF-kB", paths$Description),]$Description

nfkb.genes <- sapply(nfkb, 
                    function(x) {
                        v <- stringr::str_split(x, "/")
                                         })
nfkb.genes.list <- c(unique(unlist(nfkb.genes)))
#nfkb.genes <- nfkb.genes[-grep("MUC12|B3GNT6|ADAM|IL|CHST4", nfkb.genes)]



data.sub <- df[which(df$SYMBOL %in% nfkb.genes.list),]
melt.data <- reshape::melt(as.data.frame(data.sub), id.vars = "SYMBOL")
melt.data$variable <- gsub("\\_[0-9]*$", "", melt.data$variable)
melt.data <- tidyr::separate(melt.data, col = variable, into = c('treatment', 'hr'), sep = "-")
#melt.data <- dplyr::left_join(melt.data, nfkb.genes, by = 'SYMBOL')
#melt.data$hr <- as.numeric(melt.data$hr)


library(magrittr)
## calculate mean and stats
data.mean <- dplyr::group_by(melt.data, SYMBOL, hr) %>%
    dplyr::summarise(stdev = sd(value),
                     num = n(),
                     iqr = IQR(value),
                     min = min(value),
                     max = max(value),
		     mean = mean(value),
                     median = median(value))

data.mean$sem <- data.mean$stdev/sqrt(data.mean$num)

## claculate zscore
scale_this <- function(x) as.vector(scale(x))

melt.data$hr <- as.numeric(melt.data$hr)
data.scaled <- dplyr::group_by(melt.data, SYMBOL) %>%
    dplyr::mutate(zscore = scale_this(value)) %>%
    dplyr::group_by(SYMBOL, hr) %>%
    dplyr::summarize(mean_zscore = mean(zscore))

data.cor <- dplyr::group_by(melt.data, SYMBOL) %>%
    dplyr::summarize(cor = cor(hr, value))

top.cor <- data.cor[data.cor$cor > quantile(data.cor$cor, 0.75) | data.cor$cor < quantile(data.cor$cor, 0.25),]

## data.scaled$category <- "Glycotransferases"
## data.scaled[grep("MUC", data.scaled$SYMBOL),]$category <- "Mucins"
data.scaled <- data.scaled[which(data.scaled$SYMBOL %in% top.cor$SYMBOL),]
## data.scaled <- data.scaled[order(data.scaled$hr,-data.scaled$mean_zscore),]
## data.scaled$SYMBOL <- factor(data.scaled$SYMBOL,
##                          levels = unique(data.scaled$SYMBOL))
data.scaled$hr <- as.character(data.scaled$hr)
data.scaled$hr <- factor(data.scaled$hr, levels = c("0", "24", "48", "96"))
dat <- dplyr::select(data.scaled, hr, mean_zscore) %>% tidyr::spread(hr, mean_zscore)
ord <- hclust(dist(dat[,2:5], method = "euclidean"), method = "ward.D")$order

data.scaled$SYMBOL <- factor(data.scaled$SYMBOL, levels = unique(data.scaled$SYMBOL)[ord])

## plot
library(ggplot2)
source("ggplot2-themes.R")

figure5s1a <- ggplot(data.scaled,
              aes(y = SYMBOL, x = factor(hr))) +
    geom_tile(stat = "identity", aes(fill = mean_zscore)) +
   # facet_grid(category ~ ., scales = "free_y", space = "free", switch = "y") +
    scale_fill_distiller(name = "Z-score ", palette = "RdYlBu") +
    scale_y_discrete(position = "right") +
    ylab("") + xlab("") + 
    theme1 + 
    theme(strip.text =  element_text(size = 32),
          legend.position = "bottom",
	  legend.title = element_text(size = 32),
	  legend.key.size = unit(1.5,"cm"),
	  panel.spacing = unit(2, "lines"),
	  panel.border = element_blank()) +
	  coord_fixed(ratio =0.5) + ggtitle("A")
pdf(file = "../figures/figure5/figure5-NFkB_supplement.pdf", width = 2500/300, height = 4500/300, onefile = FALSE)
print(figure5s1a)
dev.off()

png(filename = "../figures/figure5/figure5-NFkB_supplement.png", width = 550, height = 1000)
print(figure5s1a)
dev.off()
ggsave(filename = "../figures/figure5/figure5-NFkB_supplement.eps", 
       plot = figure5s1a, 
       width = 12, height = 20)

## load data and set directory for output
## plotting

library(ggplot2)
library(ggstance)
source("ggplot2-themes.R")
library(RColorBrewer)
red.set <- brewer.pal(n = 8, name = "Reds")

## multi-volcano ---------------------------------------------------------------
data.dir <- "../results/ECOR2_hypoxia_nfkb"
library(magrittr)
ecor2.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hk.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-HK_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hypoxia.pbs <- readr::read_csv(file = file.path(data.dir,"hypoxia_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
ecor2i.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-NFKBi_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
nfkbi.pbs <- readr::read_csv(file = file.path(data.dir,"NFkBi_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hki.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-HK_over_ECOR2-HK-NFKBi.csv")) %>% dplyr::rename(SYMBOL = X1)
hypoxiai.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-HK_over_ECOR2-HK-NFKBi.csv")) %>% dplyr::rename(SYMBOL = X1)

## add comparison ID
ecor2.pbs$comparison <- "E. coli"
hk.pbs$comparison <- "HK-E. coli"
hypoxia.pbs$comparison <- "1% O2"
ecor2i.pbs$comparison <- "E. coli + SC-514"
nfkbi.pbs$comparison <- "SC-514 alone"
hki.pbs$comparison <- "HK-E. coli + SC-514"
hypoxiai.pbs$comparison <- "1% O2 + SC-514"

## bind in single dataframe
data <- rbind(ecor2.pbs,
              hk.pbs,
              hypoxia.pbs,
              ecor2i.pbs,
	      nfkbi.pbs,
	      hki.pbs,
	      hypoxiai.pbs)

plot.data <- data

plot.data$comparison <- factor(plot.data$comparison,
                               levels = rev(c("E. coli", 
                                          "1% O2",
                                          "HK-E. coli",
					  "SC-514 alone",
                                          "E. coli + SC-514",
					  "HK-E. coli + SC-514",
					  "1% O2 + SC-514")))

plot.data$status <- ifelse(plot.data$padj > 0.05 | is.na(plot.data$padj), "a",
                    ifelse(plot.data$log2FoldChange > 0, "b", "c"))

plot.data <- plot.data[order(plot.data$status),]

library(ggplot2)
source("ggplot2-themes.R")
figure5s1c <- ggplot(data = plot.data, aes(y = comparison, x = log2FoldChange)) +
    geom_point(position = position_jitter(w = 0.33), aes(fill = status, color = status), shape = 21) +
    scale_fill_manual(values = c("grey70", color.set[1], color.set[2])) +
    scale_color_manual(values = c("grey70", color.set[1], color.set[2])) +
    xlim(c(-5,5)) +
    ylab("") +
    xlab(expression(paste("log"[2],"FC over PBS"))) +
    theme1 + 
    theme(axis.text.x = element_text(size = 20)) + ggtitle("C")

print(figure5s1c)

data.dir <- "../results/ECOR2_hypoxia_nfkb"
library(magrittr)
df <- readr::read_csv(file = file.path(data.dir,"NFkBi_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)

library(clusterProfiler)
library(org.Hs.eg.db)
## SC-514-regualted genes
up.ids <- bitr(df[df$padj < 0.05 & df$log2FoldChange > 0,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db")

## background gene set
all.ids <- bitr(df$SYMBOL,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")

## over-representation test, Up-regualted
go.ec <- enrichGO(gene = up.ids$ENTREZID,
                  universe = all.ids$ENTREZID,
                  OrgDb = "org.Hs.eg.db",
                  ont = "BP",
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

## KEGG over-representation test
kegg.ec <- enrichKEGG(gene = up.ids$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)
library(ReactomePA)
## KEGG over-representation test
reactome.ec <- enrichPathway(gene = up.ids$ENTREZID,
                     # organism = 'hsa',
                      pvalueCutoff = 0.05, readable = TRUE)

## add column annotating source
reactome.ec@result$DATABASE <- "REACTOME"
kegg.ec@result$DATABASE <- "KEGG"
go.ec@result$DATABASE <- "GO"

## bind all over-representation test results
enrich.result.ec <- rbind(reactome.ec@result,
                          go.ec@result,
			  kegg.ec@result)

write.csv(enrich.result.ec, file = file.path(data.dir, "nfkbi_up_pathways.csv"))

## SC-514-regualted genes
up.ids <- bitr(df[df$padj < 0.05 & df$log2FoldChange < 0,]$SYMBOL,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db")

## background gene set
all.ids <- bitr(df$SYMBOL,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")

## over-representation test, Up-regualted
go.ec <- enrichGO(gene = up.ids$ENTREZID,
                  universe = all.ids$ENTREZID,
                  OrgDb = "org.Hs.eg.db",
                  ont = "BP",
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

## KEGG over-representation test
kegg.ec <- enrichKEGG(gene = up.ids$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)
library(ReactomePA)
## KEGG over-representation test
reactome.ec <- enrichPathway(gene = up.ids$ENTREZID,
                     # organism = 'hsa',
                      pvalueCutoff = 0.05, readable = TRUE)

## add column annotating source
reactome.ec@result$DATABASE <- "REACTOME"
kegg.ec@result$DATABASE <- "KEGG"
go.ec@result$DATABASE <- "GO"

## bind all over-representation test results
enrich.result.ec <- rbind(reactome.ec@result,
                          go.ec@result,
			  kegg.ec@result)
write.csv(enrich.result.ec, file = file.path(data.dir, "nfkbi_down_pathways.csv"))

## import data
data.dir <- "../results/ECOR2_hypoxia_nfkb"
go.up <- readr::read_csv(file = file.path(data.dir, "nfkbi_up_pathways.csv"))
go.up$direction <- "Up-regulated by SC-514"
go.down <- readr::read_csv(file = file.path(data.dir, "nfkbi_down_pathways.csv"))
go.down$direction <- "Down-regulated by SC-514"
## combine datasets
go.data <- rbind(go.up, go.down)

library(magrittr)
go.data %<>% dplyr::select(Description, qvalue, DATABASE, direction) %>%
    dplyr::group_by(DATABASE, direction) %>%
    dplyr::mutate(pctile = ecdf(-log10(qvalue))(-log10(qvalue))) %>%
    dplyr::arrange(-pctile) %>%
    subset(.$pctile > 0.9)

go.data <- go.data[order(go.data$qvalue, decreasing = TRUE),]

go.data$Description <- factor(go.data$Description,
                                levels = unique(go.data$Description))

library(ggplot2)
library(ggstance)
source("ggplot2-themes.R")
figure5s1d <- ggplot(data = go.data, aes(x = -log10(qvalue),
                                   y = Description,
                                   fill = -log10(qvalue))) +
    geom_barh(stat = "identity") +
    facet_grid(direction ~ ., scales = "free_y", space = "free") +
    theme1 +
    ylab("") +
    scale_fill_gradient(expression(paste("-log"[10],"(P-value)")),
                        high = blue.set[8],
                        low = blue.set[4]) +
    theme(strip.text.x =  element_text(size = 24, face = "bold"),
          strip.text.y =  element_text(size = 32),
          strip.background = element_rect(color = "grey30", fill = "white", size = 2.5),
          panel.background = element_rect(color = NA, fill = "grey90"),
          legend.position = "bottom",
	  legend.title = element_text(size = 24),
          axis.text.y = element_text(size = 24),
	  axis.title = element_text(size = 24),
	  plot.subtitle = element_text(size = 26, hjust = 0.5, face = "bold"),
	  panel.spacing.x = unit(0.25, "lines"),
	  panel.spacing = unit(1.5, "lines"),
	  panel.border = element_blank(),
	  legend.key.size = unit(1,"cm")) + ggtitle("D")

print(figure5s1d)           

## import datsets
library(magrittr)
gs1 <-readr::read_csv(file = "../results/supplemental/Figure5_Gene_set1.csv")$x
gs2 <-readr::read_csv(file = "../results/supplemental/Figure5_Gene_set2.csv")$x
nfkbi.dwn <- readr::read_csv(file = "../results/ECOR2_hypoxia_nfkb/NFkBi_over_PBS.csv")
nfkbi.dwn <- subset(nfkbi.dwn, nfkbi.dwn$padj < 0.05 & nfkbi.dwn$log2FoldChange < 0)$X1

library(VennDiagram)
source("ggplot2-themes.R")
venn.plot <- venn.diagram(list(gs1, gs2, nfkbi.dwn), NULL,
                          fill = c(color.set[1], color.set[2], color.set[3]),
                          alpha = 0.5,
                          lwd = 5,
                          cex = 4,
                          fontfamily = "sans",
                          lty = 0,
                          cat.font.family,
                          cat.fontfamily = "sans",
                          cat.cex = 2.5,
                          cat.pos = c(-30,30, 180),
                          cat.dist = c(0.08, 0.08,0.08),
                          category.names = c("Gene Set I\n(contact induced)",
                                              "Gene Set II\n(hypoxia-induced)",
                                              "Genes down-regulated\nby SC-514 at baseline"),
                          main.cex = 4,
                          main.fontfamily = "sans",
                          main.pos =c(0.5,0.025),
                          scaled = TRUE, euler.d = TRUE)

library(gridExtra)
library(grid)
figure5s1e <- grid.arrange(gTree(children=venn.plot), ncol = 1,
                         top = textGrob("E", hjust = 7,
                         gp = gpar(fontsize = 50, font =2)))

print(figure5s1e)

data.dir <- "../results/ECOR2_hypoxia_nfkb"
library(magrittr)
df <- readr::read_csv(file = file.path(data.dir,"NFkBi_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)

library(clusterProfiler)
library(org.Hs.eg.db)
## SC-514-regualted genes
gs1.ids <- bitr(intersect(gs1, nfkbi.dwn),
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db")

## background gene set
all.ids <- bitr(df$SYMBOL,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")

## over-representation test, Up-regualted
go.ec <- enrichGO(gene = gs1.ids$ENTREZID,
                  universe = all.ids$ENTREZID,
                  OrgDb = "org.Hs.eg.db",
                  ont = "BP",
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

## KEGG over-representation test
kegg.ec <- enrichKEGG(gene = gs1.ids$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)
library(ReactomePA)
## KEGG over-representation test
reactome.ec <- enrichPathway(gene = gs1.ids$ENTREZID,
                     # organism = 'hsa',
                      pvalueCutoff = 0.05, readable = TRUE)

## add column annotating source
reactome.ec@result$DATABASE <- "REACTOME"
kegg.ec@result$DATABASE <- "KEGG"
go.ec@result$DATABASE <- "GO"

## bind all over-representation test results
enrich.result.ec <- rbind(reactome.ec@result,
                          go.ec@result,
			  kegg.ec@result)
enrich.result.ec$gene_set <- "Gene Set I - baseline NFkB suppression"
write.csv(enrich.result.ec, file = file.path(data.dir, "nfkbi_gs1_pathways.csv"))

data.dir <- "../results/ECOR2HIO-RNAseq/"
library(magrittr)
df <- readr::read_csv(file = file.path(data.dir,"NFkBi_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)

library(clusterProfiler)
library(org.Hs.eg.db)
## SC-514-regualted genes
gs2.ids <- bitr(intersect(gs2, nfkbi.dwn),
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db")

## background gene set
all.ids <- bitr(df$SYMBOL,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")

## over-representation test, Up-regualted
go.ec <- enrichGO(gene = gs2.ids$ENTREZID,
                  universe = all.ids$ENTREZID,
                  OrgDb = "org.Hs.eg.db",
                  ont = "BP",
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

## KEGG over-representation test
kegg.ec <- enrichKEGG(gene = gs2.ids$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)
library(ReactomePA)
## KEGG over-representation test
reactome.ec <- enrichPathway(gene = gs2.ids$ENTREZID,
                     # organism = 'hsa',
                      pvalueCutoff = 0.05, readable = TRUE)

## add column annotating source
reactome.ec@result$DATABASE <- "REACTOME"
kegg.ec@result$DATABASE <- "KEGG"
go.ec@result$DATABASE <- "GO"

## bind all over-representation test results
enrich.result.ec <- rbind(reactome.ec@result,
                          go.ec@result,
			  kegg.ec@result)
enrich.result.ec$gene_set <- "Gene Set II - baseline NFkB suppression"
write.csv(enrich.result.ec, file = file.path(data.dir, "nfkbi_gs2_pathways.csv"))

## import data
data.dir <- "../results/ECOR2_hypoxia_nfkb/"
go.up <- readr::read_csv(file = file.path(data.dir, "nfkbi_gs1_pathways.csv"))
go.down <- readr::read_csv(file = file.path(data.dir, "nfkbi_gs2_pathways.csv"))
## combine datasets
go.data <- rbind(go.up, go.down)

library(magrittr)
go.data %<>% dplyr::select(Description, qvalue, DATABASE, gene_set) %>%
    dplyr::group_by(DATABASE, gene_set) %>%
    dplyr::mutate(pctile = ecdf(-log10(qvalue))(-log10(qvalue))) %>%
    dplyr::arrange(-pctile) %>%
    subset(.$pctile > 0.9)

go.data <- go.data[order(go.data$qvalue, decreasing = TRUE),]

go.data$Description <- factor(go.data$Description,
                                levels = unique(go.data$Description))

go.data$gene_set <- gsub(pattern = " - baseline NFkB suppression",
                       replacement = "",go.data$gene_set)
library(ggplot2)
library(ggstance)
source("ggplot2-themes.R")
figure5s1f <- ggplot(data = go.data, aes(x = -log10(qvalue),
                                   y = Description,
				   fill = gene_set)) +
#                                   fill = -log10(qvalue))) +
    geom_barh(stat = "identity") +
    facet_grid(gene_set ~ c("Subset supressed by\nSC-514 at baseline"), scales = "free_y", space = "free") +
    theme1 +
    ylab("") +
    scale_fill_brewer(palette = "Set1", direction = 1) +
   # scale_fill_gradient(expression(paste("-log"[10],"(P-value)")),
   #                     high = blue.set[8],
   #                     low = blue.set[4]) +
    theme(strip.text.x =  element_text(size = 24, face = "bold"),
          strip.text.y =  element_text(size = 32),
          strip.background = element_rect(color = "grey30", fill = "white", size = 2.5),
          panel.background = element_rect(color = NA, fill = "#B6D8B5"),
#          legend.position = "bottom",
	  legend.title = element_text(size = 24),
          axis.text.y = element_text(size = 24),
	  axis.title = element_text(size = 24),
	  plot.subtitle = element_text(size = 26, hjust = 0.5, face = "bold"),
	  panel.spacing.x = unit(0.25, "lines"),
	  panel.spacing = unit(1.5, "lines"),
	  panel.border = element_blank(),
	  legend.key.size = unit(1,"cm")) +
    ggtitle("F")

print(figure5s1f)            

## Figure 5 - Supplement 1 Multipanel -------------------------------------------
library(ggplot2)
library(gridExtra)
source("ggplot2-themes.R")
source("custom_fun.R")

figure5s1b <- png2ggplot("../data/figure5/p-65.png") +
    img.theme + ggtitle("B") + coord_fixed(ratio = 0.5)




layout <- rbind(c(1,2,2,4,4,4,4),
                c(1,3,3,4,4,4,4),
                c(1,5,5,5,6,6,6),
                c(1,5,5,5,6,6,6))


## PDF output
pdf(file = "../figures/supplemental/figure5s1_multipanel.pdf", width = 9500/300, height = 7500/300, onefile = FALSE)
gridExtra::grid.arrange(figure5s1a,figure5s1b,figure5s1c,figure5s1d,figure5s1e,figure5s1f,
             layout_matrix = layout)
dev.off()

## PNG output
png(filename = "../figures/supplemental/figure5s1_multipanel.png", width = 3000, height = 2500)
gridExtra::grid.arrange(figure5s1a,figure5s1b,figure5s1c,figure5s1d,figure5s1e,figure5s1f,
             layout_matrix = layout)
dev.off()

## Figure 5 - Supplement 3
## import raw data
data <- readr::read_csv(file = "../data/figure5/170825_ELISA_results.csv")

## load libraries
library(scales)
library(grid)
library(ggplot2)
source("ggplot2-themes.R")
library(magrittr)

## statistical tests
## BD-1
bd1.p <- t.test(data[data$treatment == "PBS" & data$cytokine == "BD-1",]$pg_ml,
                data[data$treatment == "E. coli" & data$cytokine == "BD-1",]$pg_ml,
                alternative = "greater")$p.value
bd1.p2 <- t.test(data[data$treatment == "PBS" & data$cytokine == "BD-1",]$pg_ml,
                data[data$treatment == "Heat-inactivated + hypoxia" & data$cytokine == "BD-1",]$pg_ml,
                alternative = "two.sided")$p.value
bd1.p3 <- t.test(data[data$treatment == "PBS" & data$cytokine == "BD-1",]$pg_ml,
                data[data$treatment == "hypoxia" & data$cytokine == "BD-1",]$pg_ml,
                alternative = "two.sided")$p.value
## BD-2
bd2.p <- t.test(data[data$treatment == "PBS" & data$cytokine == "BD-2",]$pg_ml,
                data[data$treatment == "E. coli" & data$cytokine == "BD-2",]$pg_ml,
                alternative = "less")$p.value
bd2.p2 <- t.test(data[data$treatment == "PBS" & data$cytokine == "BD-2",]$pg_ml,
                data[data$treatment == "Heat-inactivated + hypoxia" & data$cytokine == "BD-2",]$pg_ml,
                alternative = "less")$p.value
## IL-6
il6.p <- t.test(data[data$treatment == "PBS" & data$cytokine == "IL-6",]$pg_ml,
                data[data$treatment == "E. coli" & data$cytokine == "IL-6",]$pg_ml,
                alternative = "less")$p.value
il6.p2 <- t.test(data[data$treatment == "Heat-inactivated + hypoxia" & data$cytokine == "IL-6",]$pg_ml,
                data[data$treatment == "E. coli" & data$cytokine == "IL-6",]$pg_ml,
                alternative = "less")$p.value
il6.p3 <- t.test(data[data$treatment == "PBS" & data$cytokine == "IL-6",]$pg_ml,
                data[data$treatment == "Heat-inactivated + hypoxia" & data$cytokine == "IL-6",]$pg_ml,
                alternative = "less")$p.value
il6.p4 <- t.test(data[data$treatment == "PBS" & data$cytokine == "IL-6",]$pg_ml,
                data[data$treatment == "Heat-inactivated" & data$cytokine == "IL-6",]$pg_ml,
                alternative = "less")$p.value
il6.p5 <- t.test(data[data$treatment == "PBS" & data$cytokine == "IL-6",]$pg_ml,
                data[data$treatment == "hypoxia" & data$cytokine == "IL-6",]$pg_ml,
                alternative = "less")$p.value
## IL-8
il8.p <- t.test(data[data$treatment == "PBS" & data$cytokine == "IL-8",]$pg_ml,
                data[data$treatment == "E. coli" & data$cytokine == "IL-8",]$pg_ml,
                alternative = "less")$p.value
il8.p2 <- t.test(data[data$treatment == "Heat-inactivated + hypoxia" & data$cytokine == "IL-8",]$pg_ml,
                data[data$treatment == "E. coli" & data$cytokine == "IL-8",]$pg_ml,
                alternative = "less")$p.value
il8.p3 <- t.test(data[data$treatment == "PBS" & data$cytokine == "IL-8",]$pg_ml,
                data[data$treatment == "Heat-inactivated + hypoxia" & data$cytokine == "IL-8",]$pg_ml,
                alternative = "less")$p.value
il8.p4 <- t.test(data[data$treatment == "PBS" & data$cytokine == "IL-8",]$pg_ml,
                data[data$treatment == "Heat-inactivated" & data$cytokine == "IL-8",]$pg_ml,
                alternative = "less")$p.value
il8.p5 <- t.test(data[data$treatment == "PBS" & data$cytokine == "IL-8",]$pg_ml,
                data[data$treatment == "hypoxia" & data$cytokine == "IL-8",]$pg_ml,
                alternative = "less")$p.value

## VEGF
vegf.p <- t.test(data[data$treatment == "Heat-inactivated" & data$cytokine == "VEGF",]$pg_ml,
                data[data$treatment == "E. coli" & data$cytokine == "VEGF",]$pg_ml,
                alternative = "less")$p.value


data <- dplyr::group_by(data, treatment, cytokine) %>%
    dplyr::summarise(avg = mean(pg_ml),
                     num = n(),
                     sem = sd(pg_ml, na.rm = TRUE) / n(),
                     total = sum(pg_ml, na.rm = TRUE))

data$treatment <- factor(data$treatment,
                         levels = c("PBS",
                                    "E. coli",
                                    "Heat-inactivated",
                                    "hypoxia",
                                    "Heat-inactivated + hypoxia")
                         )



elisa <- ggplot(data,
                aes(x = treatment, y = avg,
                    fill = factor(treatment),
                    color = factor(treatment))) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymax = avg + sem,
                      ymin = avg - sem),
                  width = 0, color = "grey70", size = 2) +
    facet_wrap(~cytokine, scales = "free_y") +
    xlab("") +
    ylab("pg/ml") +
    scale_fill_brewer(palette = "Set1", direction = 1) +
    scale_color_brewer(palette = "Set1", direction = 1) +
    ## statistical annotations
    ## BD-1
    geom_segment(data = data[data$cytokine == "BD-1",],
                 aes(x = 1,
                     xend = 2, 
                     y = 53+4,
                     yend = 53+4),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "BD-1",],
                 aes(x = 1.5,
                     y = 55+4),
              color = "grey30",
              label = paste("P = ", format(bd1.p, digits = 1))) +
    geom_segment(data = data[data$cytokine == "BD-1",],
                 aes(x = 1,
                     xend = 5, 
                     y = 73+4,
                     yend = 73+4),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "BD-1",],
                 aes(x = 3,
                     y = 75+4),
              color = "grey30",
              label = paste("P = ", format(bd1.p2, digits = 1))) +
    geom_segment(data = data[data$cytokine == "BD-1",],
                 aes(x = 1,
                     xend = 4, 
                     y = 73,
                     yend = 73),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "BD-1",],
                 aes(x = 2.5,
                     y = 75),
              color = "grey30",
              label = paste("P = ", format(bd1.p3, digits = 1))) +
    ## BD-2
    geom_segment(data = data[data$cytokine == "BD-2",],
                 aes(x = 1,
                     xend = 2, 
                     y = 750,
                     yend = 750),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "BD-2",],
                 aes(x = 1.5,
                     y = 750+20),
              color = "grey30",
              label = paste("P = ", format(bd2.p, digits = 1))) +
    geom_segment(data = data[data$cytokine == "BD-2",],
                 aes(x = 1,
                     xend = 5, 
                     y = 750+40,
                     yend = 750+40),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "BD-2",],
                 aes(x = 3,
                     y = 750+60),
              color = "grey30",
              label = paste("P = ", format(bd2.p2, digits = 1))) +
    ## IL-6
    geom_segment(data = data[data$cytokine == "IL-6",],
                 aes(x = 1,
                     xend = 2, 
                     y = 160,
                     yend = 160),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "IL-6",],
                 aes(x = 1.5,
                     y = 160+10),
              color = "grey30",
              label = paste("P = ", format(il6.p, digits = 1))) +
    geom_segment(data = data[data$cytokine == "IL-6",],
                 aes(x = 1,
                     xend = 3, 
                     y = 160+20,
                     yend = 160+20),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "IL-6",],
                 aes(x = 2,
                     y = 160+30),
              color = "grey30",
              label = paste("P = ", format(il6.p4, digits = 1))) +
    geom_segment(data = data[data$cytokine == "IL-6",],
                 aes(x = 1,
                     xend = 4, 
                     y = 160+40,
                     yend = 160+40),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "IL-6",],
                 aes(x = 2.5,
                     y = 160+50),
              color = "grey30",
              label = paste("P = ", format(il6.p5, digits = 1))) +
    geom_segment(data = data[data$cytokine == "IL-6",],
                 aes(x = 1,
                     xend = 5, 
                     y = 160+60,
                     yend = 160+60),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "IL-6",],
                 aes(x = 3,
                     y = 160+70),
              color = "grey30",
              label = paste("P = ", format(il6.p3, digits = 1))) +
    ## IL-8
    geom_segment(data = data[data$cytokine == "IL-8",],
                 aes(x = 1,
                     xend = 2, 
                     y = 2100,
                     yend = 2100),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "IL-8",],
                 aes(x = 1.5,
                     y = 2100+100),
              color = "grey30",
              label = paste("P = ", format(il8.p, digits = 1))) +
    geom_segment(data = data[data$cytokine == "IL-8",],
                 aes(x = 1,
                     xend = 4, 
                     y = 2700,
                     yend = 2700),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "IL-8",],
                 aes(x = 2.5,
                     y = 2700+100),
              color = "grey30",
              label = paste("P = ", format(il8.p5, digits = 1))) +
    geom_segment(data = data[data$cytokine == "IL-8",],
                 aes(x = 1,
                     xend = 5, 
                     y = 3800,
                     yend = 3800),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "IL-8",],
                 aes(x = 3,
                     y = 3800+100),
              color = "grey30",
              label = paste("P = ", format(il8.p3, digits = 1))) +
    ## VEGF
    geom_segment(data = data[data$cytokine == "VEGF",],
                 aes(x = 2,
                     xend = 3, 
                     y = 650,
                     yend = 650),
                 size = 1,
		 color = "grey30")+
    geom_text(data = data[data$cytokine == "VEGF",],
                 aes(x = 2.5,
                     y = 650+25),
              color = "grey30",
              label = paste("P = ", format(vegf.p, digits = 1))) +
    theme1 +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1),
          strip.text = element_text(color = "white", size = 24),
	  strip.background = element_rect(color = "white",
                                          fill = "grey30"))
pdf(file = "../figures/supplemental/sfigure5_supp3.pdf", width = 18, height = 15, onefile = FALSE)
print(elisa)
dev.off()

## FIGURE 5 - Supplement 2 -----------------------------------------------------
## import data
data2 <- read.csv(file = "../data/figure6/170831_BD_invitro_growth/170831_OD600.csv",
                   header = TRUE, skip = 2, stringsAsFactors = FALSE)

plate2 <- read.csv(file = "../data/figure6/170831_BD_invitro_growth/170831_plate_layout.csv",
                  header = TRUE, stringsAsFactors = FALSE)

## long format
data2 <- reshape2::melt(data2, id.vars =c("Time", "Temperature..C."))
data2 <- plyr::rename(data2, c("variable"="cell"))

## reformat time to hours
data2$Time <- sapply(strsplit(data2$Time, ":"), function(x) {
    x <- as.numeric(x)
    x[1] + ((x[2] + (x[3]/60))/60)
}
)

## add cell labels from plate layout      
data2 <- plyr::join(data2, plate2, by ="cell")
data2 <- data2[complete.cases(data2),]


library(magrittr) 
plot.data <- data2 %>%
    dplyr::group_by(treatment, strain, Time) %>%
    dplyr::summarize(avg = mean(value),
                     sem = sd(value)/sqrt(n()))

library(ggplot2)
source("ggplot2-themes.R")
sfigure6.2 <- ggplot(data = plot.data[plot.data$Time < 6,], aes(x = Time, y = avg)) +
    geom_point(shape = 21, size = 5, aes(fill = factor(treatment), color = factor(treatment))) +
    geom_errorbar(aes(ymax = avg + sem, ymin = avg - sem, 
                      color = factor(treatment)), width = 0) +
    xlab("Time (h)") +
    ylab(expression(paste(Delta,"OD"[600]))) +
    scale_fill_brewer(name = "Treatment", palette = "Set1", direction = 1) +
    scale_color_brewer(name = "Treatment", palette = "Set1", direction = 1) +
    facet_wrap(~strain, ncol = 2) +
    theme1 +
    theme(legend.position = c(0.15, 0.9),
          legend.title = element_text(size = 24),
	  legend.key.size = unit(1,"cm"),
	  legend.text = element_text(size = 24),
	  strip.text = element_text(size = 32)) +
    ggtitle("A")

png(filename = "../figures/figure6/sfigure6_2.png", width = 1200, height = 800)
print(sfigure6.2)
dev.off()
ggsave(filename = "../figures/figure6/eps/sfigure6_2.eps", 
plot = sfigure6.2, 
width = 24, height = 16)

library(ggplot2)
source("ggplot2-themes.R")
figure6e <- ggplot(data = gc.data, aes(x = stringr::str_wrap(treatment, width = 20), y = K)) +
    geom_boxplot(aes(color = treatment), size = 2, outlier.size = 3) +
    facet_wrap(~strain, ncol = 2) +	
    xlab("") +
    ylab("Carrying capacity (K)") +
    scale_fill_brewer(name = "Treatment", palette = "Set1", direction = 1) +
    scale_color_brewer(name = "Treatment", palette = "Set1", direction = 1) +
    theme1 +
    theme(legend.position = "none",
          legend.title = element_text(size = 24),
	  legend.key.size = unit(1.5,"cm"),
	  legend.text = element_text(size = 24),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	  strip.text = element_text(size = 32)) +
    geom_text(data = gc.data[gc.data$strain == "E. coli ECOR2",],
                 aes(x = 1,
                     y = 0.205),
              size = 10,
              color = "grey30",
              label = paste("P = ", format(gc.test1, digits = 1))) +
    geom_text(data = gc.data[gc.data$strain == "E. coli K12",],
              aes(x = 1,
                  y = 0.220),
              size = 10,
              color = "grey30",
              label = paste("P = ", format(gc.test2, digits = 1))) +
    ggtitle("B")

png(filename = "../figures/figure6/sfigure6_2b.png", width = 1000, height = 800)
print(figure6e)
dev.off()
ggsave(filename = "../figures/figure6/eps/sfigure6_2b.png", 
plot = figure6e, 
width = 16, height = 16)

## import data 
data <- readr::read_csv(file = "../data/figure6/170825_ECOR2_BD2_invitro/9_1_2017.csv", col_names = TRUE, skip = 12)
## retrieve letter
data$letter <- substr(data$Plate_ID,1,1)

## groups dataframe
groups <- readr::read_csv(file = "../data/figure6/170825_ECOR2_BD2_invitro/170901_plate_layout.csv", col_names = TRUE)

## map to data frame
data <- dplyr::left_join(data, groups, by = 'letter')

## statistical tests
tt1 <- t.test(data[data$letter == "A" | data$letter == "D",]$Count_per_ml,
              data[data$letter == "B",]$Count_per_ml,
              alternative = 'greater')

tt2 <- t.test(data[data$letter == "A" | data$letter == "D",]$Count_per_ml,
              data[data$letter == "C",]$Count_per_ml,
              alternative = 'greater')

library(ggplot2)

source("ggplot2-themes.R")

sfig5.2c <- ggplot(data = data[data$letter != "B",], aes(x = stringr::str_wrap(treatment, width = 20), y = Count_per_ml)) +
    geom_boxplot(aes(color = treatment), size = 3) +
    ylab("CFU/mL") +
    xlab("") +
    theme1 +
    scale_color_brewer(palette = "Set1", direction = 1) +
    geom_text(x = 1, y = 7400, size = 8, color = "grey30",
             label = paste("P =" , round(tt2$p.value, digits = 3))) +
    ylim(c(5800,8100)) +
    facet_wrap(~ c(stringr::str_wrap("E. coli ECOR2 (stationary phase)", width = 20))) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1),
          strip.text = element_text(size = 24)) +
    ggtitle("C")

print(sfig5.2c)

library(gridExtra)
layout <- rbind(c(1,1,1),
                c(2,2,3))

## PDF output
pdf(file = "../figures/figure6/figure6_supplement2.pdf", width = 6400/300, height = 6400/300, onefile = FALSE)
gridExtra::grid.arrange(sfigure6.2,figure6e,sfig5.2c, layout_matrix = layout)
dev.off()

png(filename = "../figures/figure6/sfigure6_2.png", width = 1920, height = 1600)
gridExtra::grid.arrange(sfigure6.2,figure6e,sfig5.2c, layout_matrix = layout)
dev.off()

ggsave(filename = "../figures/figure6/sfigure6_2.eps", 
plot = gridExtra::grid.arrange(sfigure6.2,figure6e, sfig5.2c, layout_matrix = layout),
       width = 48, height = 12)
