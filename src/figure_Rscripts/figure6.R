## FIGURE 6 --------------------------------------------------------------------
## Figure 6A: RNA-seq mucus expression and glycotransferases timecourse
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

## FIGURE 6 --------------------------------------------------------------------
## Figure 6A: RNA-seq mucus expression and glycotransferases timecourse
## import data
## Load dataset from file
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"
library(magrittr)
df <- readr::read_csv(file = file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv")) %>% dplyr::rename(SYMBOL = X1)

## list of genes to extract
paths <- readr::read_csv("../results/ECOR2_hypoxia_nfkb/plotted-nfkb_complete-goANDreactome-results.csv")

mucins <- paths[grep("O-linked", paths$Description),]$geneID
muc.genes <- sapply(mucins, 
                    function(x) {
                        v <- stringr::str_split(x, "/")
                                         })
muc.genes <- c(unique(unlist(muc.genes))[-c(21:24)], "MUC5AC", "FUT2", "FUT3")
muc.genes <- muc.genes[-grep("MUC12|B3GNT6|ADAM|IL|CHST4", muc.genes)]



data.sub <- df[which(df$SYMBOL %in% muc.genes),]
melt.data <- reshape::melt(as.data.frame(data.sub), id.vars = "SYMBOL")
melt.data$variable <- gsub("\\_[0-9]*$", "", melt.data$variable)
melt.data <- tidyr::separate(melt.data, col = variable, into = c('treatment', 'hr'), sep = "-")
#melt.data <- dplyr::left_join(melt.data, muc.genes, by = 'SYMBOL')
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

data.scaled <- dplyr::group_by(melt.data, SYMBOL) %>%
    dplyr::mutate(zscore = scale_this(value)) %>%
    dplyr::group_by(SYMBOL, hr) %>%
    dplyr::summarize(mean_zscore = mean(zscore))

data.scaled$category <- "Glycotransferases"
data.scaled[grep("MUC", data.scaled$SYMBOL),]$category <- "Mucins"

data.scaled <- data.scaled[order(data.scaled$hr,data.scaled$mean_zscore),]
data.scaled$SYMBOL <- factor(data.scaled$SYMBOL,
                         levels = unique(data.scaled$SYMBOL))

## plot
library(ggplot2)
source("ggplot2-themes.R")

figure6a <- ggplot(data.scaled,
              aes(y = SYMBOL, x = hr)) +
    geom_tile(stat = "identity", aes(fill = mean_zscore)) +
    facet_grid(category ~ ., scales = "free_y", space = "free", switch = "y") +
    scale_fill_distiller(name = "Z-score ", palette = "RdYlBu") +
    scale_y_discrete(position = "right") +
    ylab("") + xlab("") + 
    theme1 + 
    theme(strip.text =  element_text(size = 32),
          legend.position = "bottom",
	  legend.title = element_text(size = 32),
	  legend.key.size = unit(1,"cm"),
	  panel.spacing = unit(2, "lines"),
	  panel.border = element_blank()) +
    ggtitle("A")

png(filename = "../figures/figure6/figure6a.png", width = 600, height = 1000)
print(figure6a)
dev.off()
ggsave(filename = "../figures/figure6/eps/figure6a.eps", 
       plot = figure6a, 
       width = 12, height = 20)

## FIGURE 6 --------------------------------------------------------------------
## Figure 6E: Mucin expression is NF-kB dependent
## import data
data.dir <- "../results/ECOR2_hypoxia_nfkb/"
library(magrittr)
ecor2.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hk.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-HK_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hypoxia.pbs <- readr::read_csv(file = file.path(data.dir,"hypoxia_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
ecor2i.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-NFKBi_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)

## directory where data is stored
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"

## read in the data from DESeq2 output csv files

hr48 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_48hr.csv")) %>% dplyr::rename(SYMBOL = X1)
data.dir <- "../results/ECOR2_hypoxia_nfkb/"
## add comparison ID
ecor2.pbs$comparison <- "E. coli (24h)"
hk.pbs$comparison <- "HK-E. coli"
hypoxia.pbs$comparison <- "hypoxia"
ecor2i.pbs$comparison <- "E. coli + SC-514"
hr48$comparison <- "E. coli (48h)"

## bind in single dataframe
data <- rbind(ecor2.pbs,
              hk.pbs,
              hypoxia.pbs,
              ecor2i.pbs,
	      hr48)

plot.data <- data

## list of genes to extract
paths <- readr::read_csv("../results/ECOR2_hypoxia_nfkb/plotted-nfkb_complete-goANDreactome-results.csv")

mucins <- paths[grep("O-linked", paths$Description),]$geneID
muc.genes <- sapply(mucins, 
                    function(x) {
                        v <- stringr::str_split(x, "/")
                                         })
muc.genes <- c(unique(unlist(muc.genes))[-c(21:24)], "MUC5AC", "FUT2", "FUT3")
muc.genes <- muc.genes[-grep("MUC12|B3GNT6|ADAM|IL|CHST4", muc.genes)]



plot.data <- plot.data[which(plot.data$SYMBOL %in% muc.genes),]


library(magrittr)

## claculate zscore
scale_this <- function(x) as.vector(scale(x))

data.scaled <- dplyr::group_by(plot.data, SYMBOL) %>%
    dplyr::mutate(zscore = scale_this(log2FoldChange)) %>%
    dplyr::group_by(SYMBOL, comparison) %>%
    dplyr::summarize(mean_zscore = mean(zscore))

data.scaled$category <- "Glycotransferases"
data.scaled[grep("MUC", data.scaled$SYMBOL),]$category <- "Mucins"

data.scaled$comparison <- factor(data.scaled$comparison,
                               levels = c("E. coli (48h)",
			       "E. coli (24h)", 
                                          "HK-E. coli",
                                          "E. coli + SC-514",
					  "hypoxia"))

## plot
library(ggplot2)
source("ggplot2-themes.R")

figure6e <- ggplot(data.scaled,
              aes(y = SYMBOL, x = comparison)) +
    geom_tile(stat = "identity", aes(fill = mean_zscore)) +
    facet_grid(category ~ ., scales = "free_y", space = "free", switch = "y") +
    scale_fill_distiller(name = "Z-score ", palette = "RdYlBu") +
    scale_y_discrete(position = "right") +
    ylab("") + xlab("") + 
    theme1 + 
    theme(strip.text =  element_text(size = 32),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 20),
          legend.position = "bottom",
	  legend.title = element_text(size = 32),
	  legend.key.size = unit(1,"cm"),
	  panel.spacing = unit(2, "lines"),
	  panel.border = element_blank()) +
    ggtitle("E")

png(filename = "../figures/figure6/figure6e.png", width = 600, height = 1000)
print(figure6e)
dev.off()
ggsave(filename = "../figures/figure6/eps/figure6e.eps", 
plot = figure6e, 
width = 16, height = 16)

## Figure 6 multipanel ---------------------------------------------------------
source("ggplot2-themes.R")
source("custom_fun.R")
library(ggplot2)
library(gridExtra)

## control & E. coli 10X PAS-Alcian blue
figure6b <- png2ggplot("../figures/figure6/figure6b.png") +
    img.theme + coord_fixed(ratio = 1/1.75) + ggtitle("B")

## control & E. coli 100X H&E, AB, PAS, and PAS-Alcian blue
figure6c <- png2ggplot("../figures/figure6/figure6c.png") +
    img.theme + coord_fixed(ratio = 1/2) + ggtitle("C")

## control & E. coli 60X confocal
figure6d <- png2ggplot("../figures/figure6/figure6d.png") +
    img.theme + coord_fixed(ratio = 1/2) + ggtitle("D")

figure6f <- png2ggplot("../figures/figure6/figure6f.png") +
    img.theme + coord_fixed(ratio = 1/3) + ggtitle("F")

layout <- rbind(c(1,2,2,2),
                c(3,3,4,4),
		c(5,6,6,6))

## PDF output
pdf(file = "../figures/figure6/figure6_multipanel.pdf", width = 7000/300, height = 8250/300, onefile = FALSE)
gridExtra::grid.arrange(figure6a, figure6b, figure6c, figure6d, figure6e, figure6f, layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/figure6/eps/figure6_multipanel.eps", 
       plot = gridExtra::grid.arrange(figure6a, figure6b, figure6c, figure6d, figure6e, figure6f, layout_matrix = layout), 
       width = 20, height = 20)

## PNG output
png(filename = "../figures/figure6/figure6_multipanel.png", width = 1865, height = 2000)
gridExtra::grid.arrange(figure6a, figure6b, figure6c, figure6d, figure6e, figure6f, layout_matrix = layout)
dev.off()
