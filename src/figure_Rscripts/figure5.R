## FIGURE 5 --------------------------------------------------------------------
## Figure 5B: Figure 5B: Identification of gene sets 1-4 with scatter plots
## import data
data.dir <- "../results/ECOR2_hypoxia_nfkb/"
library(magrittr)
ecor2.nfkbi <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_ECOR-NFKBi.csv")) %>% dplyr::rename(SYMBOL = X1)
hk.nfkbi <- readr::read_csv(file = file.path(data.dir,"ECOR2-HK_over_ECOR2-HK-NFKBi.csv")) %>% dplyr::rename(SYMBOL = X1)
hypoxia.nfkbi <- readr::read_csv(file = file.path(data.dir,"hypoxia_over_hypoxia-NFKBi.csv")) %>% dplyr::rename(SYMBOL = X1)

## add comparison ID
## remember the signs get reversed in the plot - labels reflect the plot orientation
ecor2.nfkbi$comparison <- "live+NFKBi vs live"#"E. coli -/+ SC-514" #"E. coli/E. coli + SC-514"
hk.nfkbi$comparison <- "dead+NFKBi vs dead"
hypoxia.nfkbi$comparison <- "hypoxia+NFKBi vs hypoxia"


df <- rbind(hk.nfkbi, hypoxia.nfkbi)


## subset columns of interest
ecor2 <- ecor2.nfkbi %>% dplyr::select(SYMBOL, log2FoldChange, padj, comparison) %>%
    dplyr::rename(ecor2_log2 = log2FoldChange, ecor2_padj = padj, comparison2 = comparison)
df <- df %>% dplyr::select(SYMBOL, log2FoldChange, padj, comparison) %>%
    dplyr::rename(hk_log2 = log2FoldChange, hk_padj = padj)

## join dataframes
df <- dplyr::left_join(df, ecor2, by = 'SYMBOL') 

## subset data prior to plotting
df1 <- subset(df, df$ecor2_padj < 0.05 |
                     df[df$comparison == "dead+NFKBi vs dead",]$hk_padj < 0.05 |
                                                                              df[df$comparison == "hypoxia+NFKBi vs hypoxia",]$hk_padj < 0.05, na.rm = TRUE)

df1$comparison <- factor(df1$comparison,
                        levels = c("dead+NFKBi vs dead",
                                   "hypoxia+NFKBi vs hypoxia"))

library(ggplot2)
source("ggplot2-themes.R")
plot.set <- df1[(df1$hk_log2 > 0 & df1$ecor2_log2 > 0 & (df1$ecor2_padj < 0.05 | df1$hk_padj < 0.05)),]

gs1 <- subset(plot.set,plot.set$comparison == "dead+NFKBi vs dead")$SYMBOL
gs2 <- subset(plot.set,plot.set$comparison == "hypoxia+NFKBi vs hypoxia")$SYMBOL

figure5b1 <- ggplot(data = df1, aes(x = -hk_log2, y = -ecor2_log2)) +
    geom_point(shape = 1, size = 3, fill = "grey", color = "grey") +
    geom_point(data = df1[(df1$hk_log2 > 0 & df1$ecor2_log2 > 0 & (df1$ecor2_padj < 0.05 | df1$hk_padj < 0.05)),],
               shape = 21, size = 3, fill = color.set[2], color = color.set[2]) + 
    stat_density_2d(data = df1[(df1$hk_log2 > 0 & df1$ecor2_log2 > 0 & (df1$ecor2_padj < 0.05 | df1$hk_padj < 0.05)),],
                    aes(fill = ..level..), geom = "polygon", na.rm = TRUE) +
    facet_grid(comparison2 ~ comparison) +
    xlim(c(-3,3)) + ylim(c(-3,3)) +
    scale_fill_distiller("Density", palette = "Spectral") +
    labs(x = expression(Log[2]*"FC"),
         y = expression(Log[2]*"FC")) +
     annotate("text",
              x = -1.5, y = -3,
              label = c(
                  paste("Gene Set I -",
                      length(
                          rownames(
                              subset(plot.set,plot.set$comparison == "dead+NFKBi vs dead"))),
                      "genes"),
                  paste("Gene Set II -",
                      length(
                          rownames(
                              subset(plot.set,plot.set$comparison == "hypoxia+NFKBi vs hypoxia"))),
                      "genes")
                  ),
              color = "grey10",
              size = 12, fontface = "bold", hjust = 0.3) +
    theme1 +
    theme(axis.title = element_text(size = 32),
          panel.spacing = unit(1, "lines"),
          #panel.border = element_rect(fill = NA, color = "white"),
          strip.text =  element_text(size = 32))

png(filename = "../figures/figure5/figure5b1.png", width = 1200, height = 600)
print(figure5b1)
dev.off()

ecor2.nfkbi <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_ECOR-NFKBi.csv")) %>% dplyr::rename(SYMBOL = X1)
hk.nfkbi <- readr::read_csv(file = file.path(data.dir,"ECOR2-HK_over_ECOR2-HK-NFKBi.csv")) %>% dplyr::rename(SYMBOL = X1)
hypoxia.nfkbi <- readr::read_csv(file = file.path(data.dir,"hypoxia_over_hypoxia-NFKBi.csv")) %>% dplyr::rename(SYMBOL = X1)

## add comparison ID
ecor2.nfkbi$comparison <- "E. coli/E. coli + SC-514"
hk.nfkbi$comparison <- "HK-E. coli/HK-E. coli + SC-514"
hypoxia.nfkbi$comparison <- "1% O2/1% O2 + SC-514"

## subset columns of interest
ecor2 <- ecor2.nfkbi %>% dplyr::select(SYMBOL, log2FoldChange, padj) %>%
    dplyr::rename(ecor2_log2 = log2FoldChange, ecor2_padj = padj)
hk <- hk.nfkbi %>% dplyr::select(SYMBOL, log2FoldChange, padj) %>%
    dplyr::rename(hk_log2 = log2FoldChange, hk_padj = padj)
hypoxia <- hypoxia.nfkbi %>% dplyr::select(SYMBOL, log2FoldChange, padj) %>%
    dplyr::rename(hypoxia_log2 = log2FoldChange, hypoxia_padj = padj)

## join dataframes
df <- dplyr::left_join(ecor2, hk, by = 'SYMBOL') %>%
    dplyr::left_join(hypoxia, by = 'SYMBOL')

## subset data prior to plotting
df <- subset(df, df$ecor2_padj < 0.05 & (df$hk_padj < 0.05 | df$hypoxia_padj < 0.05), na.rm = TRUE)
test <- df

ecor2.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hk.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-HK_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hypoxia.pbs <- readr::read_csv(file = file.path(data.dir,"hypoxia_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)

## add comparison ID
ecor2.pbs$comparison <- "E. coli/PBS"
hk.pbs$comparison <- "dead E. coli/PBS"
hypoxia.pbs$comparison <- "hypoxia/PBS"

df3 <- rbind(hk.pbs, hypoxia.pbs)

## subset columns of interest
ecor2 <- ecor2.pbs %>% dplyr::select(SYMBOL, log2FoldChange, padj, comparison) %>%
    dplyr::rename(ecor2_log2 = log2FoldChange, ecor2_padj = padj, comparison2 = comparison)
df3 <- df3 %>% dplyr::select(SYMBOL, log2FoldChange, padj, comparison) %>%
    dplyr::rename(hk_log2 = log2FoldChange, hk_padj = padj)

## join dataframes
df2 <- dplyr::left_join(df3, ecor2, by = 'SYMBOL') 


df2 <- df2[-which(df2$SYMBOL %in% df$SYMBOL),]

df <- df2

## subset data prior to plotting
df <- subset(df, df$ecor2_padj < 0.05 |
                     df[df$comparison == "dead E. coli/PBS",]$hk_padj < 0.05 |

                                                              df[df$comparison == "hypoxia/PBS",]$hk_padj < 0.05, na.rm = TRUE)

df$comparison <- factor(df$comparison,
                        levels = c("dead E. coli/PBS",
                                   "hypoxia/PBS"))

library(ggplot2)
source("ggplot2-themes.R")
#df <- df[complete.cases(df),]
plot.set <- df[(df$hk_log2 > 0 & df$ecor2_log2 > 0 & (df$ecor2_padj < 0.05 | df$hk_padj < 0.05)),]
#plot.set <- df[(df$hk_log2 > 0 & df$ecor2_log2 > 0),]

plot.set <- plot.set[complete.cases(plot.set),]

gs3 <- subset(plot.set,plot.set$comparison == "dead E. coli/PBS")$SYMBOL
gs4 <- subset(plot.set,plot.set$comparison == "hypoxia/PBS")$SYMBOL

figure5b2 <- ggplot(data = df, aes(x = hk_log2, y = ecor2_log2)) +
    geom_point(shape = 1, size = 3, fill = "grey", color = "grey") +
     ## geom_point(data = df[(df$hk_log2 < 0 & df$ecor2_log2 < 0 ),],
     ##            shape = 21, size = 3, fill = color.set[2], color = color.set[2]) + 
     ## stat_density_2d(data = df[(df$hk_log2 < 0 & df$ecor2_log2 < 0),],
     ##                 aes(fill = ..level..), geom = "polygon", na.rm = TRUE) +
     geom_point(data = plot.set,
                shape = 21, size = 3, fill = color.set[2], color = color.set[2]) + 
     stat_density_2d(data = plot.set,
                     aes(fill = ..level..), geom = "polygon", na.rm = TRUE) +

     facet_grid(comparison2 ~ comparison) +
    xlim(c(-3,3)) + ylim(c(-3,3)) +
    scale_fill_distiller("Density", palette = "Spectral") +
    labs(x = expression(Log[2]*"FC"),
         y = expression(Log[2]*"FC")) +
    annotate("text",
             x = -1.5, y = 3,
              label = c(
                  paste("Gene Set III -",
                      length(
                          rownames(
                              subset(plot.set,plot.set$comparison == "dead E. coli/PBS"))),
                      "genes"),
                  paste("Gene Set IV -",
                      length(
                          rownames(
                              subset(plot.set,plot.set$comparison == "hypoxia/PBS"))),
                      "genes")
                  ),

              color = "grey10",
              size = 12, fontface = "bold", hjust = 0.1) +
    theme1 +
    theme(axis.title = element_text(size = 32),
          panel.spacing = unit(1, "lines"),
          #panel.border = element_rect(fill = NA, color = "white"),
          strip.text =  element_text(size = 32))

png(filename = "../figures/figure5/figure5b2.png", width = 1200, height = 600)
print(figure5b2)
dev.off()

library(gridExtra)
png(filename = "../figures/figure5/figure5b.png", width = 2400, height = 600)
gridExtra::grid.arrange(figure5b1, figure5b2, ncol = 2)
dev.off()


ggsave(filename = "../figures/figure5/eps/figure5b.eps", 
       plot = gridExtra::grid.arrange(figure5b1, figure5b2, ncol = 2), 
       width = 48, height = 12)

## Output supplemental table with genes in each set
write.csv(gs1, file = "../results/supplemental/Figure5_Gene_set1.csv")
write.csv(gs2, file = "../results/supplemental/Figure5_Gene_set2.csv")
write.csv(gs3, file = "../results/supplemental/Figure5_Gene_set3.csv")
write.csv(gs4, file = "../results/supplemental/Figure5_Gene_set4.csv")

## FIGURE 5 --------------------------------------------------------------------
## Figure 5C: GO & REACTOME plot
## import data
## load data
data.dir <- "../results/ECOR2_hypoxia_nfkb/"
data <- readr::read_csv(file = file.path(data.dir, "nfkb-dependent_complete-goANDreactome-results.csv.bak"))

data2 <- readr::read_csv(file = file.path(data.dir, "nfkb-independent_complete-goANDreactome-results.csv"))

reactome.data <- readr::read_csv(file = file.path(data.dir, "REACTOME-requireslive_nfkbindependent_ecoli_subset.csv"))

reactome.data$comparison <- "live-hypoxia"

data2 <- rbind(data2, reactome.data)
data$nfkb <- "SC-514 suppressed"
data2$nfkb <- ""
data2$comparison <- gsub(", nfkb-independent", "", data2$comparison)
data <- rbind(data, data2)

data$comparison <- paste0(data$comparison,"\n",data$nfkb)
## load list of pathways to plot
select <- readr::read_csv(file = file.path(data.dir, "unique_nfkb_pathways.csv"), col_names = TRUE)

## subset data to selections
#data <- data[which(data$Description %in% unique(select$X2)),]
data <- data[which(data$Description %in% select$Description),]
test <- data #check
data <- subset(data, data$pvalue < 0.004 | data$Description == "Defensins")
data <- dplyr::left_join(data, select, by = "Description")

data$Description <- gsub("Hypoxia-inducible Factor Alpha", "HIF1a",data$Description)
data$Description <- gsub("kappa", "k",data$Description)
data$Description <- gsub("toll-like receptor", "TLR",data$Description)
data$Description <- gsub("Toll Like Receptor", "TLR",data$Description)
data$Description <- gsub("Toll-Like Receptors", "All TLR",data$Description)
data$Description <- gsub(" initiated on plasma membrane", "",data$Description)
data$Description <- gsub("interleukin", "IL",data$Description)
data$Description <- gsub("positive regulation", "up-regulation",data$Description)
data$Description <- gsub("Mediated ", "",data$Description)
#data$Description <- gsub("response to ", "",data$Description)
data$Description <- gsub("mediated ", "",data$Description)
data$Description <- gsub("negative regulation", "down-regulation",data$Description)
data$Description <- gsub("nucleotide-binding oligomerization domain containing", "NOD-like",data$Description)
data$Description <- gsub("lic process", "",data$Description)
data$Description <- gsub("catabo", "catabolism",data$Description)
data$Description <- gsub("metabo", "metabolism",data$Description)
data$Description <- gsub("signaling pathway", "pathway",data$Description)
data$Description <- gsub("TLR T", "T",data$Description)
data$Description <- gsub("The citric acid ", "",data$Description)
data$Description <- gsub("\\(TCA\\)", "TCA",data$Description)
data$Description <- gsub("TAK1 activates NFkB by phosphorylation and activation of IKKs complex", "TAK1 activates NFkB",data$Description)
data$Description <- gsub("\\(TLR5\\)|\\(TLR10\\)|\\(TLR9\\)|\\(TLR2\\)|\\(TLR3\\)|\\(TLR4\\)|", "",data$Description)
data$Description <- gsub("  ", " ",data$Description)

data$comparison <- gsub(pattern = "live-HK",
                       replacement = "Associated with\nbacterial contact", data$comparison)
data$comparison <- gsub(pattern = "live-hypoxia",
                       replacement = "Associated with\nbacterial hypoxia", data$comparison)
data$category <- gsub("Innate and adaptive defense", "Innate and adaptive\ndefense",data$category)

## apply gene ratio cut-off
data <- data[data$gene_ratio > 0.01,]

data$comparison <- factor(data$comparison,
                          levels = c(
                              "Associated with\nbacterial contact\nSC-514 suppressed",
                              "Associated with\nbacterial hypoxia\nSC-514 suppressed",
                              "Associated with\nbacterial contact\n",                 
                              "Associated with\nbacterial hypoxia\n",                 
                              "hk-hypoxia\nSC-514 suppressed",                        
                              "hk-hypoxia\n"))   

 data$category <- factor(data$category,
                         levels = c(
                             "TLR signaling",
                             "Innate and adaptive\ndefense",
                             "NF-kB signaling",
                             "Epithelial barrier",
                             "Mucins",
                             "Metabolism",
                             "Mitochondria",
                             "Hypoxia",
                             "Development",
                             "Cell cycle"))

data <- data[order( data$comparison,data$gene_ratio, decreasing = FALSE),]

data$Description <- factor(data$Description,
                                levels = unique(data$Description))


#data$combo_group <- paste0(data$hypoxia_status," \n ",data$hk_status)
write.csv(data, file = file.path(data.dir, "plotted-nfkb_complete-goANDreactome-results.csv"))

library(ggplot2)
library(ggstance)
source("ggplot2-themes.R")

## Use "Gene Set" names
data$comparison <- gsub("Associated with\nbacterial contact\nSC-514 suppressed", "Gene Set I", data$comparison)
data$comparison <- gsub("Associated with\nbacterial hypoxia\nSC-514 suppressed", "Gene Set II", data$comparison)
data$comparison <- gsub("Associated with\nbacterial contact\n", "Gene Set III", data$comparison)
data$comparison <- gsub("Associated with\nbacterial hypoxia\n", "Gene Set IV", data$comparison)                 


figure5c <- ggplot(data = data[data$comparison != "hk-hypoxia\n" &
                                        #   data$category != "Hypoxia" &
                                             data$category != "NF-kB signaling" &                               
                                              data$comparison != "hk-hypoxia\nSC-514 suppressed",],
                   aes(x = gene_ratio*100, y = Description, fill = -log10(pvalue))) +
    geom_barh(stat = "identity") +
    facet_grid(category ~ comparison, scales = "free_y", space = "free_y") +
    scale_fill_gradient(expression(paste("-log"[10],"(P-value)")),
                        high = blue.set[8],
                        low = blue.set[4])+
    xlab("% genes matched to pathway") + ylab("") +
    theme1 + 
    theme(strip.text.x =  element_text(size = 32, face = "bold"),
          strip.text.y =  element_text(size = 32, 
                                       angle = 360),
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
	  legend.key.size = unit(1,"cm")) + ggtitle("C")


png(filename = "../figures/figure5/figure5c.png", width = 2400, height = 1550)
print(figure5c)
dev.off()
ggsave(filename = "../figures/figure5/eps/figure5c.eps", 
       plot = figure5c, 
       width = 48, height = 31)

## multipanel plot -------------------------------------------------------------
library(ggplot2)
source("ggplot2-themes.R")
source("custom_fun.R")

if (file.exists("../figures/figure5/figure5a-cairo.Rdata") == TRUE) {
    load(file = "../figures/figure5/figure5a-cairo.Rdata")
} else {    
    ##https://www.stat.auckland.ac.nz/~paul/Reports/Rlogo/Rlogo.html
    ##convertPicture("../figures/figure5/figure5a.svg", "../figures/figure5/figure5a-cairo.svg")
    ## this step takes a while
    figure5a <- grImport2::readPicture("../figures/figure5/figure5a-cairo.svg")
    save(figure5a, file = "../figures/figure5/figure5a-cairo.Rdata")
    load(file = "../figures/figure5/figure5a-cairo.Rdata")
}

library(grid)
library(gridSVG)
library(grImport2)
fig4a <- gTree(children = gList(pictureGrob(figure5a, ext = "gridSVG")))

figure5a <- qplot(1:100, 1:100, alpha = I(0)) +
    theme_bw() +
    annotation_custom(fig4a, xmin = -Inf,
                      xmax = Inf,
                      ymin = -Inf,
                      ymax = Inf) +
    img.theme + ggtitle("A") + coord_fixed(ratio = 0.354)

layout <- rbind(c(1),
                c(1),
                c(2),
		c(3),
                c(3),
                c(3),
                c(3))

library(gridExtra)
library(grid)

figure5b <- gridExtra::grid.arrange(figure5b1, figure5b2, ncol = 2,
                                  top = textGrob("B", hjust = 32,
                                                 gp = gpar(fontsize = 45, font = 2)))

## PDF output
pdf(file = "../figures/figure5/figure5_multipanel.pdf", width = 9200/300, height = 14562/300, onefile = FALSE)
gridExtra::grid.arrange(figure5a, figure5b, figure5c, layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/figure5/eps/figure5_multipanel.eps", 
       plot = gridExtra::grid.arrange(figure5a, figure5b, figure5c, layout_matrix = layout), 
       width = 48, height = 67, limitsize = FALSE, device = "eps")

## PNG output
png(filename = "../figures/figure5/figure5_multipanel.png", width = 2400, height = 3350)
gridExtra::grid.arrange(figure5a, figure5b, figure5c, layout_matrix = layout)
dev.off()
