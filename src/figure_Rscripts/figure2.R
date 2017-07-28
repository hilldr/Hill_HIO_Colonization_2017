## FIGURE 2 --------------------------------------------------------------------
## Figure 2A: Multi-volcano plot
## import data
## directory where data is stored
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"

## read in the data from DESeq2 output csv files
library(magrittr)
hr24 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_24hr.csv")) %>% dplyr::rename(SYMBOL = X1)
hr48 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_48hr.csv")) %>% dplyr::rename(SYMBOL = X1)
hr96 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_96hr.csv")) %>% dplyr::rename(SYMBOL = X1)

## add hr variable
hr24$hr <- 24
hr48$hr <- 48
hr96$hr <- 96

## bind in single dataframe
data <- rbind(hr24,
              hr48,
              hr96)

## make a copy to edit for plotting
plot.data <- data
## create status catergory for assigning colors
plot.data$status <- ifelse(plot.data$padj > 0.05 | is.na(plot.data$padj), "a",
                    ifelse(plot.data$log2FoldChange > 0, "b", "c"))
## sort by status
plot.data <- plot.data[order(plot.data$status),]

## generate plot
library(ggplot2)
source("ggplot2-themes.R")
figure2a <- ggplot(data = plot.data, aes(x = factor(hr), y = log2FoldChange)) +
    geom_point(position = position_jitter(w = 0.33), aes(fill = status, color = status), shape = 21) +
    scale_fill_manual(values = c("grey70", color.set[1], color.set[2])) +
    scale_color_manual(values = c("grey70", color.set[1], color.set[2])) +
    ylim(c(-5,5)) +
    xlab("Hours post-microinjection") +
    ylab(expression(paste("-log"[2],"(HIO + E. coli / HIO + PBS)"))) +
    theme1 +
    ggtitle("A")

dir.create(path = "../figures/figure2/eps", recursive = TRUE)

png(filename = "../figures/figure2/figure2a.png", width = 800, height = 1000)
print(figure2a)
dev.off()

ggsave(filename = "../figures/figure2/eps/figure2a.eps", 
       plot = figure2a, 
       width = 16, height = 20)

## FIGURE 2 --------------------------------------------------------------------
## Figure 2B: PCA
## import data
## read in table with sample metadata
samples <- readr::read_csv(file = "../data/RNA-seq/sample_key.csv")

## subsetting rules
## exclude ES cells (for now)
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

## Load dataset from file
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"
df <- readr::read_csv(file = file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv"))

## subset to numeric columns only
num.data <- df[,sapply(df,is.numeric)]
## use colnames as group names, remove '_##"
group <- sub("_.*$", "", colnames(num.data))

## calculate variance by row (gene)
var <- apply(num.data, 1, sd, na.rm = TRUE)
## adjust cut off according to variance percentile
pca.data <- num.data[var > quantile(var, 0.9) & var != 0,]
pca <- prcomp(t(pca.data),scale = TRUE,center = TRUE)
scores <- data.frame(colnames(pca.data), pca$x[,1:ncol(pca$x)],group)
scores$short_name <- colnames(num.data)
scores <- dplyr::left_join(scores, samples, by = 'short_name')
scores$hr <- c(48,48,48,96,96,96,96,24,24,24,24,24,0,0,0,0)

## PCA plot
library(ggplot2)
library(RColorBrewer)
source("ggplot2-themes.R")
source("custom_fun.R")

figure2b <- qplot(x = PC1, y = PC2, data = scores) +  
    scale_fill_brewer(name = "hr", palette = "Spectral", direction = -1) +
    scale_color_brewer(name = "hr", palette = "Spectral", direction = -1) +  
    geom_point(shape = 21, aes(fill = factor(hr), color = factor(hr)), size = 12, stroke = 3) +
    theme1 +
    theme(legend.position="right",
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(color = "white",fill = "white"),
	  panel.border = element_rect(fill = NA, color = "grey70"),
          legend.title = element_text(size = 32),
          legend.text = element_text(size = 18)) +
    coord_fixed(ratio = 1) +
    xlab(paste("PC1 (",percent(round(summary(pca)$importance[2,1],4)),")",sep = "")) +
    ylab(paste("PC2 (",percent(round(summary(pca)$importance[2,2],4)),")",sep = "")) +
    geom_hline(yintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    ggtitle ("B")


png(filename = "../figures/figure2/figure2b.png", width = 1000, height = 1000)
print(figure2b)
dev.off()
ggsave(filename = "../figures/figure2/eps/figure2b.eps", 
       plot = figure2b, 
       width = 20, height = 20)

## FIGURE 2 --------------------------------------------------------------------
## Figure 2C: GSEA timecourse heatmap
## import data
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"

## load GO GSEA stats
hr24.gsea <- readr::read_csv(file = file.path(data.dir,"hr24.GO-GSEA.csv"))
hr48.gsea <- readr::read_csv(file = file.path(data.dir,"hr48.GO-GSEA.csv"))
hr96.gsea <- readr::read_csv(file = file.path(data.dir,"hr96.GO-GSEA.csv"))

## add column annotating comparison origin
hr24.gsea$comparison <- "24"
hr48.gsea$comparison <- "48"
hr96.gsea$comparison <- "96"
## bind all gsea and export
all.gsea1 <- rbind(hr24.gsea,
                  hr48.gsea,
                  hr96.gsea)

## load REACTOME GSEA stats
hr24.gsea <- readr::read_csv(file = file.path(data.dir,"hr24.REACTOME-GSEA.csv"))
hr48.gsea <- readr::read_csv(file = file.path(data.dir,"hr48.REACTOME-GSEA.csv"))
hr96.gsea <- readr::read_csv(file = file.path(data.dir,"hr96.REACTOME-GSEA.csv"))

## add column annotating comparison origin
hr24.gsea$comparison <- "24"
hr48.gsea$comparison <- "48"
hr96.gsea$comparison <- "96"
## bind all gsea and export
all.gsea2 <- rbind(hr24.gsea,
                  hr48.gsea,
                  hr96.gsea)
## add database
all.gsea1$Database <- "GO"
all.gsea2$Database <- "REACTOME"

all.gsea <- rbind(all.gsea1, all.gsea2)

all.gsea$new_Description <- paste0(all.gsea$Description, " (", all.gsea$Database,")")

cats <- readr::read_csv(file = "../figures/figure2/GSEA-pathways-plotted.csv")

new.gsea <- dplyr::inner_join(cats, all.gsea, by = "Description")

new.gsea <- new.gsea[order(new.gsea$Category, new.gsea$NES),]

## shorter names
new.gsea$Description <- gsub("Hypoxia-inducible Factor Alpha", "HIF1a",new.gsea$Description)
new.gsea$Description <- gsub("kappa", "k",new.gsea$Description)
new.gsea$Description <- gsub("interleukin", "IL",new.gsea$Description)
new.gsea$Description <- gsub("positive regulation", "up-regulation",new.gsea$Description)
new.gsea$Description <- gsub("Mediated ", "",new.gsea$Description)
new.gsea$Description <- gsub("mediated ", "",new.gsea$Description)
new.gsea$Description <- gsub("negative regulation", "down-regulation",new.gsea$Description)

new.gsea$new_Description <- factor(new.gsea$new_Description,
                                levels = unique(new.gsea$new_Description))

new.gsea$Description <- factor(new.gsea$Description,
                                levels = unique(new.gsea$Description))

library(ggplot2)
source("ggplot2-themes.R")

figure2c <- ggplot(new.gsea,
              aes(y = Description, x = comparison)) +
    geom_tile(stat = "identity", aes(fill = NES)) +
    facet_grid(Category ~ ., scales = "free_y", space = "free", switch = "y") +
    scale_fill_distiller(name = "NES ", palette = "RdYlBu") +
    scale_y_discrete(position = "right") +
    ylab("") + xlab("") + 
    theme1 + 
    theme(strip.text =  element_text(size = 23),
          legend.position = "bottom",
	  legend.title = element_text(size = 32),
	  legend.key.size = unit(1,"cm"),
	  panel.spacing = unit(2, "lines"),
	  panel.border = element_blank()) +
    ggtitle("C")


png(filename = "../figures/figure2/figure2c.png", width = 900, height = 1400)
print(figure2c)
dev.off()
ggsave(filename = "../figures/figure2/eps/figure2c.eps", 
       plot = figure2c, 
       width = 18, height = 28)

## FIGURE 2 --------------------------------------------------------------------
## Figure 2D: RNA expression time-course
## import data
## directory where data is stored
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"

## read in the data from DESeq2 output csv files
library(magrittr)
hr24 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_24hr.csv")) %>% dplyr::rename(SYMBOL = X1)
hr48 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_48hr.csv")) %>% dplyr::rename(SYMBOL = X1)
hr96 <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS_96hr.csv")) %>% dplyr::rename(SYMBOL = X1)

## add hr variable
hr24$hr <- 24
hr48$hr <- 48
hr96$hr <- 96

## bind in single dataframe
data <- rbind(hr24,
              hr48,
              hr96)
## genes to retrieve
gene.list <- c("DEFB1","DEFB4A","IL6", "CXCL8", "VEGFA")

## create subset
plot.data <- data[which(data$SYMBOL %in% gene.list),]

library(ggplot2)
source("ggplot2-themes.R")
plot.data$SYMBOL <- factor(plot.data$SYMBOL,
                           levels = c("DEFB1",
                                      "DEFB4A",
                                      "IL6",
                                      "CXCL8",
                                      "VEGFA"))

figure2d <- ggplot(data = plot.data,
                   aes(x = factor(hr), y = log2FoldChange, fill = log2FoldChange)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ SYMBOL, scales = "free_y", ncol = 1, strip.position = "left") +
    geom_errorbar(aes(ymax = log2FoldChange + lfcSE, ymin = log2FoldChange - lfcSE),
                  width = 0, color = "grey70", size = 1) +
    scale_fill_gradient2(name = expression(paste(Log[2],"FC")), low = color.set[2], high = color.set[1]) +

    ylab("") + xlab("hr") +
    theme1 + 
    ggtitle("D") +
    theme(strip.placement = "outside",
          strip.text = element_text(size = 24),
          legend.position = "bottom",
          legend.key.size = unit(1,"cm"),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 18))

png(filename = "../figures/figure2/figure2d.png", width = 300, height = 1000)
print(figure2d)
dev.off()
ggsave(filename = "../figures/figure2/eps/figure2d.eps", 
       plot = figure2d, 
       width = 12, height = 20)

## FIGURE 2 --------------------------------------------------------------------
## Figure 2E: Time-dependent protein secretion
## import data
## import data
data <- readr::read_csv(file = "../data/figure1/161206_survival/survival_and_ELISA.csv")
## create unique IDs
data$ID <- paste(data$well,data$plate, sep = "P")
## load libraries
library(ggplot2)
library(reshape2)
library(dplyr)

## reshape dataframe
data <- dplyr::select(data, ID, day, bd1, bd2, il6, il8, vegf, treatment, Survival = dead)
data <- dplyr::rename(data, BD1 = bd1, BD2 = bd2, IL6 = il6, IL8 = il8, VEGF = vegf)
data <- reshape2::melt(data, id.vars = c('ID', 'day', 'treatment'))

## dotplot +/- SEM
data2 <- dplyr::group_by(data,day, treatment, variable) %>%
    dplyr::summarise(avg = mean(value,na.rm = TRUE),
              sem = sd(value, na.rm = TRUE)/n(), 
              total = sum(value, na.rm = TRUE))

data2[data2$variable == "Survival",]$avg <- 1 - data2[data2$variable == "Survival",]$avg
data2[data2$variable == "Survival",]$total <- 48 - data2[data2$variable == "Survival",]$total

library(ggplot2)
source("ggplot2-themes.R")
## adjust order for legend
data2$treatment <- factor(data2$treatment,
                     levels = c("PBS", "E. coli"))


figure2e  <- ggplot(data = data2[data2$variable != "Survival",],
                    aes(x = day, y = avg, color = treatment, fill = treatment)) +
    geom_errorbar(data = data2[data2$variable != "Survival",],aes(x = day,
                      ymin = avg - sem, ymax = avg + sem, 
                      color = treatment),
                  width = 0, size = 1) +
    geom_line(aes(color = treatment), size = 0.5) +
    geom_point(data = data2[data2$variable != "Survival",], shape = 21, size = 5) +
    facet_wrap(~ variable, scales = "free_y", ncol = 1, 
               strip.position = "left") +
    xlab("day") +
    ylab("") +
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9)) +
    scale_colour_brewer(palette = "Set1", direction = -1) +
    scale_fill_brewer(palette = "Set1", direction = -1) +
    guides(fill = guide_legend(title = NULL),
           color = guide_legend(title = NULL)) +
    theme1 +
    theme(text = element_text(size = 24),
          axis.text.y = element_text(size = 16),
          legend.position = "bottom",
          strip.placement = "outside",
	  strip.text = element_text(color = "white", size = 24),
	  strip.background = element_rect(color = "white", fill = "grey30"),
          legend.key.size = unit(2,"cm"),
	  legend.text = element_text(size = 32)) +
    ggtitle("E")

png(filename = "../figures/figure2/figure2e.png", width = 600, height = 1000)
print(figure2e)
dev.off()
ggsave(filename = "../figures/figure2/eps/figure2e.eps", 
       plot = figure2e, 
       width = 12, height = 20)

## multipanel plot -------------------------------------------------------------

library(gridExtra)
layout <- rbind(c(1,1,1,2,2,2,2),
                c(3,3,3,3,4,5,5),
                c(3,3,3,3,4,5,5))                

grid.arrange(figure2a, figure2b, figure2c, figure2d, figure2e,
             layout_matrix = layout)

## PDF output
pdf(file = "../figures/figure2/figure2_multipanel.pdf", width = 6000/300, height = 7500/300, onefile = FALSE)
gridExtra::grid.arrange(figure2a, figure2b, figure2c, figure2d, figure2e, layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/figure2/eps/figure2_multipanel.eps", 
       plot = gridExtra::grid.arrange(figure2a, figure2b, figure2c, figure2d, figure2e, layout_matrix = layout), 
       width = 32, height = 40)

## PNG output
png(filename = "../figures/figure2/figure2_multipanel.png", width = 1600, height = 2000)
gridExtra::grid.arrange(figure2a, figure2b, figure2c, figure2d, figure2e, layout_matrix = layout)
dev.off()

## plot
library(ggplot2)
source("ggplot2-themes.R")

figure2s1 <- ggplot(data.scaled,
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
	  coord_fixed(ratio =0.5)
pdf(file = "../figures/figure2/figure2-NFkB_supplement.pdf", width = 2500/300, height = 4500/300, onefile = FALSE)
print(figure2s1)
dev.off()

png(filename = "../figures/figure2/figure2-NFkB_supplement.png", width = 550, height = 1000)
print(figure2s1)
dev.off()
ggsave(filename = "../figures/figure2/figure2-NFkB_supplement.eps", 
       plot = figure2s1, 
       width = 12, height = 20)
