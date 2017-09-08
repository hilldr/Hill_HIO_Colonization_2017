## FIGURE 7 --------------------------------------------------------------------
## Figure 7A: RNA-seq mucus expression and glycotransferases timecourse
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

figure7a <- ggplot(data.scaled,
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

png(filename = "../figures/figure7/figure7a.png", width = 600, height = 1000)
print(figure7a)
dev.off()
ggsave(filename = "../figures/figure7/eps/figure7a.eps", 
       plot = figure7a, 
       width = 12, height = 20)

## FIGURE 7 --------------------------------------------------------------------
## Figure 7E: Mucin expression is NF-kB dependent
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

figure7e <- ggplot(data.scaled,
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

png(filename = "../figures/figure7/figure7e.png", width = 600, height = 1000)
print(figure7e)
dev.off()
ggsave(filename = "../figures/figure7/eps/figure7e.eps", 
plot = figure7e, 
width = 16, height = 16)

## Figure 7 multipanel ---------------------------------------------------------
source("ggplot2-themes.R")
source("custom_fun.R")
library(ggplot2)
library(gridExtra)

## control & E. coli 10X PAS-Alcian blue
figure7b <- png2ggplot("../figures/figure7/figure7b.png") +
    img.theme + coord_fixed(ratio = 1/1.75) + ggtitle("B")

## control & E. coli 100X H&E, AB, PAS, and PAS-Alcian blue
figure7c <- png2ggplot("../figures/figure7/figure7c.png") +
    img.theme + coord_fixed(ratio = 1/2) + ggtitle("C")

## control & E. coli 60X confocal
figure7d <- png2ggplot("../figures/figure7/figure7d.png") +
    img.theme + coord_fixed(ratio = 1/2) + ggtitle("D")

figure7f <- png2ggplot("../figures/figure7/figure7f.png") +
    img.theme + coord_fixed(ratio = 1/3) + ggtitle("F")

layout <- rbind(c(1,2,2,2),
                c(3,3,4,4),
		c(5,6,6,6))

## PDF output
pdf(file = "../figures/figure7/figure7_multipanel.pdf", width = 7000/300, height = 8250/300, onefile = FALSE)
gridExtra::grid.arrange(figure7a, figure7b, figure7c, figure7d, figure7e, figure7f, layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/figure7/eps/figure7_multipanel.eps", 
       plot = gridExtra::grid.arrange(figure7a, figure7b, figure7c, figure7d, figure7e, figure7f, layout_matrix = layout), 
       width = 20, height = 20)

## PNG output
png(filename = "../figures/figure7/figure7_multipanel.png", width = 1865, height = 2000)
gridExtra::grid.arrange(figure7a, figure7b, figure7c, figure7d, figure7e, figure7f, layout_matrix = layout)
dev.off()
