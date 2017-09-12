## FIGURE 3 --------------------------------------------------------------------
## Figure 3A: EdU and Ki67 count data
## execute counting script
system("./figure3_count.sh")
## import data
counts <- readr::read_delim(file = "../data/figure3/cell_counts.tsv", delim = "\t",
                            col_names = FALSE, skip = 1)

## separate count data from file name
library(magrittr) # enables pipe operation ('%>%')
counts$count <- stringr::str_split_fixed(string = counts$X1, pattern = " ../data/figure3/d", n = 2)[,1] %>% as.numeric()
counts$file <- stringr::str_split_fixed(string = counts$X1, pattern = " ../data/figure3/d", n = 2)[,2]

## remove unnecessary info from files names, store as metadata string
counts$meta <- gsub("api_ki67_edu_counts/", "", counts$file) %>%
    gsub("_Counts.csv", "", .) %>%
    gsub("_\\(w\\)ecad", "", .) %>%
    gsub("\\(g\\)", "", .) %>%
    gsub("\\(r\\)", "", .) %>%
    gsub("\\(b\\)", "", .) %>%
    gsub("h", "", .) %>%
    ## correct typos in file names
    gsub("edui", "edu", .) %>%
    gsub("ki67i", "ki67", .)

## parse metadata string into data columns
## replicate number
counts$rep <- stringr::str_split_fixed(string = counts$meta, pattern = "_", n = 4)[,1] %>% as.numeric()
## treatment group
counts$treatment <- stringr::str_split_fixed(string = counts$meta, pattern = "_", n = 4)[,2]
## time
counts$hr <- stringr::str_split_fixed(string = counts$meta, pattern = "_", n = 4)[,3] %>% as.numeric()
## marker
counts$marker <- stringr::str_split_fixed(string = counts$meta, pattern = "_", n = 4)[,4]
## some instances of "edu' have extra space
counts$marker <- gsub("edu ", "edu", counts$marker)

## reformatting data (wide)
data <- counts %>% dplyr::select(treatment, hr, rep, marker, count) %>%
    tidyr::spread(key = marker, value = count)

## calculate proportions
data$edu_prop <- data$edu/data$dapi
data$edu_prop <- ifelse(data$edu_prop > 1, 1, data$edu_prop)
data$ki67_prop <- data$ki67/data$dapi

library(ggplot2)
source("ggplot2-themes.R")
library(magrittr)

data.cts <- data %>% dplyr::select(treatment, hr, edu_prop, ki67_prop) %>%
    tidyr::gather(marker, prop, 3:4)
data.cts$marker <- ifelse(data.cts$marker == "edu_prop", "EdU", "Ki67")
data.cts$treatment <- data.cts$treatment %>%
    gsub("pbs", "PBS", .) %>%
    gsub("ecoli", "E. coli", .)

test1 <- t.test(data.cts[data.cts$treatment == "PBS" & data.cts$hr == 24 & data.cts$marker == "EdU",]$prop,
                data.cts[data.cts$treatment == "E. coli" & data.cts$hr == 24 & data.cts$marker == "EdU",]$prop[-2],
                alternative = "two.sided")
test2 <- t.test(data.cts[data.cts$treatment == "PBS" & data.cts$hr == 48 & data.cts$marker == "EdU",]$prop[-4],
                data.cts[data.cts$treatment == "E. coli" & data.cts$hr == 48 & data.cts$marker == "EdU",]$prop,
                alternative = "two.sided")
test3 <- t.test(data.cts[data.cts$treatment == "PBS" & data.cts$hr == 24 & data.cts$marker == "Ki67",]$prop,
                data.cts[data.cts$treatment == "E. coli" & data.cts$hr == 24 & data.cts$marker == "Ki67",]$prop,
                alternative = "two.sided")
test4 <- t.test(data.cts[data.cts$treatment == "PBS" & data.cts$hr == 48 & data.cts$marker == "Ki67",]$prop,
                data.cts[data.cts$treatment == "E. coli" & data.cts$hr == 48 & data.cts$marker == "Ki67",]$prop,
                alternative = "two.sided")

## summary data for plotting
data.sum <- data %>% dplyr::group_by(treatment, hr) %>%
    dplyr::summarise(mean_edu = mean(edu_prop),
                     mean_ki67 = mean(ki67_prop)) %>%
    tidyr::gather(marker, mean, c('mean_edu', 'mean_ki67'))
data.sum2 <- data %>% dplyr::group_by(treatment, hr) %>%
    dplyr::summarise(sem_edu = sd(edu_prop, na.rm = TRUE)/n(),
                     sem_ki67 = sd(ki67_prop, na.rm = TRUE)/n()) %>%
    tidyr::gather(marker, sem, c('sem_edu', 'sem_ki67'))
data.sum$marker %<>% gsub("mean_edu", "EdU", .) %>%
    gsub("mean_ki67", "Ki67", .)
data.sum2$marker %<>% gsub("sem_edu", "EdU", .) %>%
    gsub("sem_ki67", "Ki67", .)
data.sum %<>% dplyr::left_join(data.sum2, by = c('marker', 'hr', 'treatment'))
data.sum$treatment %<>% gsub("pbs", "PBS", .) %>%
    gsub("ecoli", "E. coli", .) %>%
    factor(levels = c("PBS", "E. coli"))

figure3a <- ggplot(data = data.sum, aes(x = hr, y = mean)) +
    geom_smooth(aes(color = treatment)) +
    geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem, color = treatment),
                   width = 0, size = 2) +
     geom_point(shape = 21, size = 8, stroke = 2, aes(fill = treatment)) +
     facet_grid(marker~., scales = "free_y") +
     ylab("Proportion of E-cadherin(+) cells") +
     xlab("Time post-microinjection (h)") +
     scale_fill_brewer(name = "", palette = "Set1", direction = -1) +
     scale_color_brewer(name = "", palette = "Set1", direction = -1) +
     scale_x_continuous(breaks = c(0, 24, 48, 72, 96)) +
     theme1 +
     geom_text(data = data.sum[data.sum$marker == "EdU",],
                 aes(x = 26,
                     y = 0.72),
              size = 6, 
              color = "grey30",
              label = paste("P = ", round(test1$p.value, 3))) +
     geom_text(data = data.sum[data.sum$marker == "EdU",],
               aes(x = 48,
                   y = 0.22),
               size = 6, 
               color = "grey30",
               label = paste("P = ", round(test2$p.value, 3))) +
     geom_text(data = data.sum[data.sum$marker == "Ki67",],
               aes(x = 48,
                   y = 0.13),
               size = 6, 
               color = "grey30",
               label = paste("P = ", round(test4$p.value, 3))) +
     theme(strip.text = element_text(size = 36),
           legend.position = "bottom",
           legend.key.size = unit(1.5,"cm"),
           legend.title = element_text(size = 24),
           legend.text = element_text(size = 32)) +
     ggtitle("B")

png(filename = "../figures/figure3/figure3a.png", width = 800, height = 800)
print(figure3a)
dev.off()
ggsave(filename = "../figures/figure3/eps/figure3a.eps", 
plot = figure3a, 
width = 16, height = 16)

## FIGURE 3 --------------------------------------------------------------------
## Figure 3E: Cell proliferation heatmap
## load dataset from file
data.dir <- "../results/ECOR2HIO_24-96-RNAseq/"
library(magrittr)
df <- readr::read_csv(file = file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv")) %>%
    dplyr::rename(SYMBOL = X1)
## list of genes to extract
gene.list <- data.frame(SYMBOL = c("SPDEF", "MUC2", "KLF5","CCND1","MKI67", "PCNA", "SOX9", "DPP4", "LYZ", "DEFA5", "LGR5", "OLFM4", "SI", "TREH", "FABP2", "NGN3", "CHGA", "DEFA6", "HOPX", "CDH1", "LRIG1", "BMI1", "MGAM", "LCT"),
                        category = c("Goblet", "Goblet", "Stem/Progenitor","Proliferation","Proliferation", "Proliferation", "Stem/Progenitor", "Enterocyte", "Paneth", "Paneth", "Stem/Progenitor", "Stem/Progenitor", "Enterocyte", "Enterocyte", "Enterocyte", "Enteroendocrine", "Enteroendocrine", "Paneth","Stem/Progenitor", "Enterocyte", "Stem/Progenitor", "Stem/Progenitor","Enterocyte","Enterocyte"))
                        

## extract genes
data.sub <- df[which(df$SYMBOL %in% gene.list$SYMBOL),]
## melt dataframe to wide format
melt.data <- reshape::melt(as.data.frame(data.sub), id.vars = "SYMBOL")
## remove replicate numbers from column headers
melt.data$variable <- gsub("\\_[0-9]*$", "", melt.data$variable)
## extract timepoint data from column headers
melt.data <- tidyr::separate(melt.data, col = variable, into = c('treatment', 'hr'), sep = "-")
## calculate mean 
data.mean <- dplyr::group_by(melt.data, SYMBOL, hr) %>%
    dplyr::summarise(stdev = sd(value),
                     num = n(),
		     mean = mean(value),
                     median = median(value)) %>%
    dplyr::left_join(gene.list, by = 'SYMBOL')

## claculate zscore
scale_this <- function(x) as.vector(scale(x))

data.scaled <- dplyr::group_by(melt.data, SYMBOL) %>%
    dplyr::mutate(zscore = scale_this(value)) %>%
    dplyr::group_by(SYMBOL, hr) %>%
    dplyr::summarize(mean_zscore = mean(zscore)) %>%
    dplyr::left_join(gene.list, by = 'SYMBOL')

## plot
library(ggplot2)
source("ggplot2-themes.R")

## hierarchical clustering
dat <- dplyr::select(data.scaled, hr, mean_zscore) %>% tidyr::spread(hr, mean_zscore)
ord <- hclust(dist(dat[,2:5], method = "euclidean"), method = "ward.D")$order

data.scaled$SYMBOL <- factor(data.scaled$SYMBOL, levels = unique(data.scaled$SYMBOL)[ord])

figure3e <- ggplot(data.scaled,
              aes(y = SYMBOL, x = hr)) +
    geom_tile(stat = "identity", aes(fill = mean_zscore)) +
    facet_grid(category ~ ., scales = "free_y", space = "free", switch = "y") +
    scale_fill_distiller(name = "Z-score ", palette = "RdYlBu", direction = -1) +
    scale_y_discrete(position = "right") +
    ylab("") + xlab("") + 
    theme1 + 
    theme(strip.text.y =  element_text(size = 36,
                                       angle = 180, hjust = 1),
          strip.background = element_rect(color = "white", fill = "white", size = 1),
          legend.position = "bottom",
	  legend.title = element_text(size = 32),
	  legend.key.size = unit(2,"cm"),
	  panel.spacing = unit(2, "lines"),
	  panel.border = element_blank()) +
	  
    coord_fixed(ratio = 2) +
    ggtitle("E")

png(filename = "../figures/figure3/figure3e.png", width = 500, height = 1000)
print(figure3e)
dev.off()
ggsave(filename = "../figures/figure3/eps/figure3e.eps", 
plot = figure3e, 
width = 16, height = 16)

## Figure 6 multipanel ---------------------------------------------------------
source("ggplot2-themes.R")
source("custom_fun.R")
library(ggplot2)
library(gridExtra)

## EdU staining
figure3b <- png2ggplot("../data/figure3/EdU_figure_cropped.png") +
    img.theme + coord_fixed(ratio = 1/1.4) + ggtitle("A")

## Sox9 staining
figure3c <- png2ggplot("../data/figure3/sox9_figure.png") +
    img.theme + coord_fixed(ratio = 1/1.75) + ggtitle("C")

## dppIV staining
figure3d <- png2ggplot("../data/figure3/dppiv_figure.png") +
    img.theme + coord_fixed(ratio = 1/1.75) + ggtitle("D")

figure3a <- figure3a + coord_fixed(ratio = 0.5)

## layout <- rbind(c(1,1,1,2,2,2),
##                 c(1,1,1,2,2,2),
##                 c(3,3,3,3,3),
##                 c(3,3,3,3,3),
##                 c(4,4,4,4,4),
## 		c(4,4,4,4,4))

layout <- rbind(c(1,1,1,2),
                c(3,3,4,4))

## PDF output
pdf(file = "../figures/figure3/figure3_multipanel.pdf", width = 8000/300, height = 7500/300, onefile = FALSE)
gridExtra::grid.arrange(figure3b, figure3a, figure3c, figure3d, layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/figure3/eps/figure3_multipanel.eps", 
       plot = gridExtra::grid.arrange(figure3b, figure3a, figure3c, figure3d, layout_matrix = layout), 
       width = 20, height = 20)

## PNG output
png(filename = "../figures/figure3/figure3_multipanel.png", width = 1865, height = 2000)
gridExtra::grid.arrange(figure3b, figure3a, figure3c, figure3d, layout_matrix = layout)
dev.off()
