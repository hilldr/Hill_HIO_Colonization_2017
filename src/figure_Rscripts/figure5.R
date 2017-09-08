## FIGURE 5 --------------------------------------------------------------------
## Figure 5A: title
## import data
## load data
data.dir <- "../results/ECOR2_hypoxia_nfkb/"
library(magrittr)
ecor2.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
## load list of defensins
def <- readr::read_delim(file = "../figures/figure5/defensin_gene_family.txt", delim = "\t")

amps <- c("CAMP", "CAP18","FALL39", "LTF", "HAMP", "HTN1", "HTN3", "LEAP2", "PPBP", "CXCL14", "REG3A", "CCL20", "CXCL1", "CXCL2", "CXCL3", "CXCL9","CXCL10","CXCL11", "CXCL12","CXCL13", "CCL1", "CCL8", "CCL11","CCL13","CCL17","CCL18","CCL19","CCL25","CCL21", "SLPI", "CXCL4","CXCL6", "CCL27", "LYZ", "CCL24", "CCL26")
amps <- c(unlist(def[,2]), amps)
plot.data <- ecor2.pbs[which(ecor2.pbs$SYMBOL %in% amps),]
#plot.data <- plot.data[complete.cases(plot.data),]

## sort by fold change
plot.data <- plot.data[order(plot.data$log2FoldChange),]
plot.data$SYMBOL <- factor(plot.data$SYMBOL,
                         levels = unique(plot.data$SYMBOL))

## plot
library(ggplot2)
library(ggstance)
source("ggplot2-themes.R")

figure5a <- ggplot(plot.data[plot.data$baseMean > 1,],
              aes(x = log2FoldChange, y = SYMBOL, fill = log2FoldChange, color =log2FoldChange)) +
    geom_barh(stat = "identity") +
    geom_errorbarh(aes(xmax = log2FoldChange + lfcSE, xmin = log2FoldChange - lfcSE), width = 0, color = "grey70") +
    scale_fill_gradient2(low = color.set[2], high = color.set[1]) +
    scale_color_gradient2(low = color.set[2], high = color.set[1]) +
    xlab(expression(paste(Log[2],"FC")))+ ylab("") +
    geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "solid") +
    theme1 + 
    ggtitle("A", subtitle = "Antimicrobial peptides") +
    theme(plot.subtitle = element_text(size = 32, hjust = 0.5, face = "bold"))

png(filename = "../figures/figure5/figure5a.png", width = 600, height = 1000)
print(figure5a)
dev.off()
ggsave(filename = "../figures/figure5/eps/figure5a.eps", 
       plot = figure5a, 
       width = 12, height = 16)

## FIGURE 5 --------------------------------------------------------------------
## Figure 5B: DEFB4A expression is NF-kB independent
## import data
data.dir <- "../results/ECOR2_hypoxia_nfkb/"
library(magrittr)
ecor2.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hk.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-HK_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hypoxia.pbs <- readr::read_csv(file = file.path(data.dir,"hypoxia_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
ecor2i.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-NFKBi_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)

## add comparison ID
ecor2.pbs$comparison <- "E. coli"
hk.pbs$comparison <- "heat-inactivated"
hypoxia.pbs$comparison <- "hypoxia"
ecor2i.pbs$comparison <- "E. coli + SC-514"

## bind in single dataframe
data <- rbind(ecor2.pbs,
              hk.pbs,
              hypoxia.pbs,
              ecor2i.pbs)

plot.data <- data

plot.data$comparison <- factor(plot.data$comparison,
                               levels = c(
                                   "E. coli", 
                                   "heat-inactivated",
                                   "E. coli + SC-514",
                                   "hypoxia"))

data.sub <- plot.data[which(plot.data$SYMBOL %in% c("DEFB4A", "DEFB4B")),]

## t.test
source("custom_fun.R")
tt1 <- t.test2(m1 = data.sub[data.sub$comparison == 'E. coli',]$log2FoldChange[1],
               m2 = data.sub[data.sub$comparison == 'heat-inactivated',]$log2FoldChange[1],
               s1 = data.sub[data.sub$comparison == 'E. coli',]$lfcSE[1],
               s2 = data.sub[data.sub$comparison == 'heat-inactivated',]$lfcSE[1],
               n1 = 4, n2 = 4) 

tt2 <- t.test2(m1 = data.sub[data.sub$comparison == 'E. coli',]$log2FoldChange[1],
               m2 = data.sub[data.sub$comparison == 'hypoxia',]$log2FoldChange[1],
               s1 = data.sub[data.sub$comparison == 'E. coli',]$lfcSE[1],
               s2 = data.sub[data.sub$comparison == 'hypoxia',]$lfcSE[1],
               n1 = 4, n2 = 4)

tt3 <- t.test2(m1 = data.sub[data.sub$comparison == 'E. coli',]$log2FoldChange[2],
               m2 = data.sub[data.sub$comparison == 'heat-inactivated',]$log2FoldChange[2],
               s1 = data.sub[data.sub$comparison == 'E. coli',]$lfcSE[2],
               s2 = data.sub[data.sub$comparison == 'heat-inactivated',]$lfcSE[2],
               n1 = 4, n2 = 4) 

tt4 <- t.test2(m1 = data.sub[data.sub$comparison == 'E. coli',]$log2FoldChange[2],
               m2 = data.sub[data.sub$comparison == 'hypoxia',]$log2FoldChange[2],
               s1 = data.sub[data.sub$comparison == 'E. coli',]$lfcSE[2],
               s2 = data.sub[data.sub$comparison == 'hypoxia',]$lfcSE[2],
               n1 = 4, n2 = 4)

## subset to show DEFB4A only (possibly plot duplicate in supplement)
figure5b <- ggplot(data.sub[data.sub$baseMean > 1 & data.sub$SYMBOL == 'DEFB4A',],
                   aes(y = log2FoldChange, x = comparison,
                       fill = comparison, color = comparison)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymax = log2FoldChange + lfcSE, ymin = log2FoldChange - lfcSE), size = 1.5, width = 0, color = "grey70") +
    scale_fill_brewer(palette = "Spectral") +
    scale_color_brewer(palette = "Spectral") +
    facet_grid(. ~ SYMBOL) +
    ylab(expression(paste(Log[2],"FC")))+ xlab("") +
    geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "solid") +
    theme1 + 
    ggtitle("B") +
    annotate("segment",
             x = 1, xend = 2, y = 9.5, yend = 9.5,
             color = "grey30", size = 2) +  
    annotate("text", x = 1.5, y = 9.75, label = paste("P =", round(tt1[4], digits =3)), size = 10, color = "grey30") +
    theme(plot.subtitle = element_text(size = 32, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	  strip.text.x =  element_text(size = 32, face = "bold"))

png(filename = "../figures/figure5/figure5b.png", width = 600, height = 1000)
print(figure5b)
dev.off()
ggsave(filename = "../figures/figure5/eps/figure5b.eps", 
       plot = figure5b, 
       width = 12, height = 16)

## subset to show DEFB4A only (possibly plot duplicate in supplement)
defb4b.supp <- ggplot(data.sub[data.sub$baseMean > 1 & data.sub$SYMBOL == 'DEFB4B',],
                   aes(y = log2FoldChange, x = comparison,
                       fill = comparison, color = comparison)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymax = log2FoldChange + lfcSE, ymin = log2FoldChange - lfcSE), size = 1.5, width = 0, color = "grey70") +
    scale_fill_brewer(palette = "Spectral") +
    scale_color_brewer(palette = "Spectral") +
    facet_grid(. ~ SYMBOL) +
    ylab(expression(paste(Log[2],"FC")))+ xlab("") +
    geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "solid") +
    theme1 + 
    annotate("segment",
             x = 1, xend = 2, y = 6.5, yend = 6.5,
             color = "grey30", size = 2) +  
    annotate("text", x = 1.5, y = 6.75, label = paste("P =", round(tt3[4], digits =3)), size = 10, color = "grey30") +
    theme(plot.subtitle = element_text(size = 32, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	  strip.text.x =  element_text(size = 32, face = "bold"))

png(filename = "../figures/supplemental/DEFB4B-expression.png", width = 600, height = 1000)
print(defb4b.supp)
dev.off()

pdf(file = "../figures/supplemental/DEFB4B-expression.pdf", width = 2400/300, height = 4000/300)
print(defb4b.supp)
dev.off()

## FIGURE 5 --------------------------------------------------------------------
## Figure 5C: /E. coli/ colonization is associated with increased h\beta{}D-2 expression
## import data

data <- readr::read_csv(file = "../data/figure1/161206_survival/survival_and_ELISA.csv")
## create unique IDs
data$ID <- paste(data$well,data$plate, sep = "P")
## load libraries
library(reshape2)
library(dplyr)

## reshape dataframe
data <- dplyr::select(data, ID, day, bd1, bd2, il6, il8, vegf, treatment)
data <- dplyr::rename(data, BD1 = bd1, BD2 = bd2, IL6 = il6, IL8 = il8, VEGF = vegf)
data <- reshape2::melt(data, id.vars = c('ID', 'day', 'treatment'))

## convert from pg/ml to ug/ml
data$value <- (data$value)/1000 # pg/1ml = 2*pg/0.5ml
data$treatment <- gsub("HK- E. coli", "heat-inactivated", data$treatment)
data$treatment <- factor(data$treatment,
                         levels = c("PBS", "E. coli", "heat-inactivated", "E. coli + SC-514", "hypoxia", "E. coli + vehicle"),
                         ordered = TRUE)
dat <- data[data$variable == "BD2" & data$day == 1,]
stats <- t.test(dat[dat$treatment == "E. coli",]$value,dat[dat$treatment == "PBS",]$value, alternative = "two.sided")
hbd2.stats2 <- t.test(dat[dat$treatment == "E. coli",]$value,dat[dat$treatment == "E. coli + SC-514",]$value, alternative = "two.sided")
hbd2.stats3 <- t.test(dat[dat$treatment == "E. coli",]$value,dat[dat$treatment == "heat-inactivated",]$value, alternative = "two.sided")
hbd2.stats4 <- t.test(dat[dat$treatment == "E. coli",]$value,dat[dat$treatment == "hypoxia",]$value, alternative = "two.sided")

library(scales)
library(grid)
library(ggplot2)
source("ggplot2-themes.R")
data <- subset(data, data$treatment != "E. coli + vehicle")

## mean data
library(magrittr)
data_mean <- data %>% dplyr::group_by(day, treatment, variable) %>%
    summarize(mean = mean(value),
              sem = sd(value)/sqrt(n()))

figure5c <- ggplot(data[data$variable == "BD2" & data$day == 1,],
               aes(x = treatment, y = value, color = factor(treatment))) +
    geom_boxplot(size = 2, outlier.size = 3) +
    xlab("") +
    ylab(expression(paste("BD-2 (",mu,"g/mL)"))) +
    ylim(c(0, 2)) +
    scale_color_brewer(palette = "Set1") +
#    scale_fill_manual(values = c(color.set[2], color.set[1])) +
    annotate("segment", x = 1, xend = 2, y = 1.6, yend = 1.6, size = 2, color = "grey30") +
    annotate("text", x = 1.5, y = 1.65,
             label = paste("P = ",format(stats$p.value, digits = 1)),
             size = 6, color = "grey30") +
    annotate("segment", x = 2, xend = 3, y = 1.7, yend = 1.7, size = 2, color = "grey30") +
    annotate("text", x = 2.5, y = 1.75,
             label = paste("P = ",format(hbd2.stats3$p.value, digits = 1)),
             size = 6, color = "grey30") +
    annotate("segment", x = 2, xend = 4, y = 1.8, yend = 1.8, size = 2, color = "grey30") +
    annotate("text", x = 3, y = 1.85,
             label = paste("P = ",format(hbd2.stats2$p.value, digits = 1)),
             size = 6, color = "grey30") +
    annotate("segment", x = 2, xend = 5, y = 1.9, yend = 1.9, size = 2, color = "grey30") +
    annotate("text", x = 3.5, y = 1.95,
             label = paste("P = ",format(hbd2.stats4$p.value, digits = 1)),
             size = 6, color = "grey30") +
    theme1 +
    ggtitle("C") +
        theme(plot.subtitle = element_text(size = 32, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	  strip.text.x =  element_text(size = 32, face = "bold"))

figure5c.alt <- ggplot(data[data$variable == "BD2" & data$day == 1,],
               aes(x = treatment, y = value, fill = factor(treatment), color = factor(treatment))) +
    geom_bar(stat = "identity", size = 0.5) +
    geom_errorbar(aes(ymax = mean + sem, ymin = mean - sem), size = 1.5, width = 0, color = "grey70") +
    xlab("") +
    ylab(expression(paste("BD-2 (",mu,"g/mL)"))) +
    ylim(c(0, 1)) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
#    scale_fill_manual(values = c(color.set[2], color.set[1])) +
    annotate("segment", x = 1, xend = 2, y = 0.75, yend = 0.75, size = 2, color = "grey30") +
    annotate("text", x = 1.5, y = 0.78,
             label = paste("P = ",format(stats$p.value, digits = 1)),
             size = 6, color = "grey30") +
    annotate("segment", x = 2, xend = 3, y = 0.82, yend = 0.82, size = 2, color = "grey30") +
    annotate("text", x = 2.5, y = 0.85,
             label = paste("P = ",format(hbd2.stats3$p.value, digits = 1)),
             size = 6, color = "grey30") +
    annotate("segment", x = 2, xend = 4, y = 0.89, yend = 0.89, size = 2, color = "grey30") +
    annotate("text", x = 3, y = 0.92,
             label = paste("P = ",format(hbd2.stats4$p.value, digits = 1)),
             size = 6, color = "grey30") +
    annotate("segment", x = 2, xend = 5, y = 0.96, yend = 0.96, size = 2, color = "grey30") +
    annotate("text", x = 3.5, y = 0.99,
             label = paste("P = ",format(hbd2.stats2$p.value, digits = 1)),
             size = 6, color = "grey30") +
    theme1 +
    ggtitle("C") +
        theme(plot.subtitle = element_text(size = 32, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	  strip.text.x =  element_text(size = 32, face = "bold"))


png(filename = "../figures/figure5/figure5c.png", width = 600, height = 800)
print(figure5c)
dev.off()

ggsave(filename = "../figures/figure5/eps/figure5c.eps", 
       plot = figure5c, 
       width = 12, height = 16)

## FIGURE 5 --------------------------------------------------------------------
## Figure 5D: BD-2 suppresses growth of /E. coli/ /in vitro/
## import data
data2 <- read.csv(file = "../data/figure5/160518_ECOR2_BD2/160517_ECOR2_BD2_OD600.csv",
                   header = TRUE, skip = 2, stringsAsFactors = FALSE)

plate2 <- read.csv(file = "../data/figure5/160518_ECOR2_BD2/160505_OD600-ECOR2-BD1_plate.csv",
                  header = TRUE, stringsAsFactors = FALSE)
data2 <- reshape2::melt(data2, id.vars =c("Time", "Temperature..C."))
data2 <- plyr::rename(data2, c("variable"="cell"))

data2$Time <- sapply(strsplit(data2$Time, ":"), function(x) {
    x <- as.numeric(x)
    x[1] + ((x[2] + (x[3]/60))/60)
}
)
      
data2 <- plyr::join(data2, plate2, by ="cell")
data2$dose <- as.numeric(data2$dose)
data2 <- data2[complete.cases(data2),]

test1 <- t.test(data2[data2$Time == 18 & data2$dose == 1,]$value,
                data2[data2$Time == 18 & data2$dose == 0.1,]$value,
                alternative = "two.sided")

test2 <- t.test(data2[data2$Time == 18 & data2$dose == 1,]$value,
                data2[data2$Time == 18 & data2$dose == 1e-8,]$value,
                alternative = "two.sided")

plot.data <- data2[data2$dose == 1e-8| data2$dose == 0.1 |data2$dose == 1,] %>%
    dplyr::group_by(dose, Time) %>%
    dplyr::summarize(avg = mean(value), sem = sd(value)/sqrt(n()))

plot.data[plot.data$dose == 1e-8,]$dose <- 0

library(ggplot2)
source("ggplot2-themes.R")
figure5d <- ggplot(data = plot.data, aes(x = Time, y = avg)) +
    geom_point(shape = 21, size = 5, aes(fill = factor(dose), color = factor(dose))) +
    geom_errorbar(aes(ymax = avg + sem, ymin = avg - sem, 
                      color = factor(dose)), width = 0) +
    xlab("Time (h)") +
    ylab(expression(paste(Delta,"OD"[600]))) +
    scale_fill_manual(name = expression(paste("BD-2 (",mu,"g/mL)")), values = c(color.set[3], color.set[2], color.set[1])) +
    scale_color_manual(name = expression(paste("BD-2 (",mu,"g/mL)")), values = c(color.set[3], color.set[2], color.set[1])) +
    scale_x_continuous(limits = c(0,21),  breaks = c(0, 3, 6, 9, 12, 15, 18)) +
    theme1 +
    annotate("segment", x = 18.5, xend = 18.5,
             y = mean(plot.data[plot.data$Time == 18 & plot.data$dose == 0.1,]$avg),
             yend = mean(plot.data[plot.data$Time == 18 & plot.data$dose == 1,]$avg),
             size = 2, color = "grey30") +
    annotate("text", x = 20, y = 0.3,
             label = paste("P = ", round(test1$p.value, digits = 3)),
             size = 10, color = "grey30") +
    annotate("segment", x = 18.25, xend = 18.25,
             y = mean(plot.data[plot.data$Time == 18 & plot.data$dose == 0,]$avg),
             yend = mean(plot.data[plot.data$Time == 18 & plot.data$dose == 1,]$avg),
             size = 2, color = "grey30") +
    annotate("text", x = 20, y = 0.45,
             label = paste("P = ", round(test2$p.value, digits = 3)),
             size = 10, color = "grey30") +
    theme(legend.position = c(0.8, 0.15),
          legend.title = element_text(size = 24),
	  legend.key.size = unit(1.5,"cm"),
	  legend.text = element_text(size = 24)) +

    ggtitle("D")

png(filename = "../figures/figure5/figure5d.png", width = 1200, height = 800)
print(figure5d)
dev.off()
ggsave(filename = "../figures/figure5/eps/figure5d.eps", 
plot = figure5d, 
width = 24, height = 16)

## calculate growth curves & carrying capacity
library(magrittr)
## limit dataset to log - stationary phase (exclude post-stationary phase)
data3 <- subset(data2, data2$Time < 10)

gc.data <- dplyr::group_by(data3, cell) %>%
    dplyr::summarise(K = growthcurver::SummarizeGrowth(Time, value)[1]$vals$k,
                     DT = growthcurver::SummarizeGrowth(Time, value)[1]$vals$t_gen,
                     t_mid = growthcurver::SummarizeGrowth(Time, value)[1]$vals$t_mid) %>%
    dplyr::right_join(dplyr::select(data2, cell, dose) , by = "cell") %>%
    dplyr::distinct(cell, K, DT, dose, t_mid)

gc.data <- subset(gc.data, gc.data$dose > 1e-02 | gc.data$dose == 1e-8)
gc.data[gc.data$dose == 1e-8,]$dose <- 0

## statistical tests
gc.test1 <- t.test(gc.data[gc.data$dose == 0,]$K,
                   gc.data[gc.data$dose == 1,]$K,
                   alternative = "two.sided")
gc.test2 <- t.test(gc.data[gc.data$dose == 0,]$K,
                   gc.data[gc.data$dose == 0.1,]$K,
                   alternative = "two.sided")
gc.test3 <- t.test(gc.data[gc.data$dose == 0,]$DT,
                   gc.data[gc.data$dose == 1,]$DT,
                   alternative = "two.sided")
gc.test4 <- t.test(gc.data[gc.data$dose == 0,]$DT,
                   gc.data[gc.data$dose == 0.1,]$DT,
                   alternative = "two.sided")

library(ggplot2)
source("ggplot2-themes.R")
figure5e <- ggplot(data = gc.data, aes(x = factor(dose), y = K)) +
    geom_boxplot(aes(color = factor(dose)), size = 2, outlier.size = 3) +
    xlab(expression(paste("BD-2 (",mu,"g/mL)"))) +
    ylab("Carrying capacity (K)") +
    scale_fill_manual(name = expression(paste("BD-2 (",mu,"g/mL)")), values = c(color.set[3], color.set[2], color.set[1])) +
    scale_color_manual(name = expression(paste("BD-2 (",mu,"g/mL)")), values = c(color.set[3], color.set[2], color.set[1])) +
    theme1 +
    theme(legend.position = "none",
          legend.title = element_text(size = 24),
	  legend.key.size = unit(1.5,"cm"),
	  legend.text = element_text(size = 24)) +
    annotate("text", x = 2, y = 0.49,
             label = paste("P = ", round(gc.test2$p.value, digits = 3)),
             size = 6, color = "grey30") +
    annotate("text", x = 3, y = 0.47,
             label = paste("P = ", round(gc.test1$p.value, digits = 4)),
             size = 6, color = "grey30") +
    ggtitle("E")

figure5f <- ggplot(data = gc.data, aes(x = factor(dose), y = DT)) +
    geom_boxplot(aes(color = factor(dose)), size = 2, outlier.size = 3) +
    xlab(expression(paste("BD-2 (",mu,"g/mL)"))) +
    ylab("Doubling time (h)") +
    scale_fill_manual(name = expression(paste("BD-2 (",mu,"g/mL)")), values = c(color.set[3], color.set[2], color.set[1])) +
    scale_color_manual(name = expression(paste("BD-2 (",mu,"g/mL)")), values = c(color.set[3], color.set[2], color.set[1])) +
    theme1 +
    theme(legend.position = "none",
          legend.title = element_text(size = 24),
	  legend.key.size = unit(1.5,"cm"),
	  legend.text = element_text(size = 24)) +
    ggtitle("F")

png(filename = "../figures/figure5/figure5e.png", width = 800, height = 800)
#print(figure5e)
gridExtra::grid.arrange(figure5e, figure5f, ncol = 2)
dev.off()
ggsave(filename = "../figures/figure5/eps/figure5e.eps", 
plot = figure5e, 
width = 16, height = 16)

## Figure 5 multipanel ---------------------------------------------------------
source("ggplot2-themes.R")
library(ggplot2)
library(gridExtra)
source("custom_fun.R")

layout <- rbind(c(1,2,3),
                c(4,4,5))


## PDF output
pdf(file = "../figures/figure5/figure5_multipanel.pdf", width = 6000/300, height = 6000/300, onefile = FALSE)
gridExtra::grid.arrange(figure5a, figure5b, figure5c, figure5d, figure5e, layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/figure5/eps/figure5_multipanel.eps", 
       plot = gridExtra::grid.arrange(figure5a, figure5b, figure5c, figure5d, figure5e, layout_matrix = layout), 
       width = 20, height = 20)

## PNG output
png(filename = "../figures/figure5/figure5_multipanel.png", width = 1600, height = 1600)
gridExtra::grid.arrange(figure5a, figure5b, figure5c, figure5d, figure5e, layout_matrix = layout)
dev.off()

## calculate growth curves & carrying capacity
library(magrittr)
## limit dataset to log - stationary phase (exclude post-stationary phase)
data3 <- subset(data2, data2$Time < 6)

gc.data <- dplyr::group_by(data3, cell) %>%
    dplyr::summarise(K = growthcurver::SummarizeGrowth(Time, value)[1]$vals$k,
                     DT = growthcurver::SummarizeGrowth(Time, value)[1]$vals$t_gen,
                     t_mid = growthcurver::SummarizeGrowth(Time, value)[1]$vals$t_mid) %>%
    dplyr::right_join(dplyr::select(data2, cell, treatment, strain) , by = "cell") %>%
    dplyr::distinct(cell, K, DT, treatment, strain, t_mid)

## statistical tests
gc.test1 <- t.test(gc.data[gc.data$treatment == 'PBS' & gc.data$strain == "E. coli ECOR2",]$K,
                   gc.data[gc.data$treatment == '1.0 ug/ml BD-2' & gc.data$strain == "E. coli ECOR2",]$K,
                   alternative = "greater")$p.value
gc.test2 <- t.test(gc.data[gc.data$treatment == 'PBS' & gc.data$strain == "E. coli K12",]$K,
                   gc.data[gc.data$treatment == '1.0 ug/ml BD-2' & gc.data$strain == "E. coli K12",]$K,
                   alternative = "greater")$p.value
