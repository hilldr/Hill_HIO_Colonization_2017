## FIGURE 8 --------------------------------------------------------------------
## Figure 8A: Regulation of cytoskeleton organization and cell-cell adhesion is NF-kB dependent
## import data
data.dir <- "../results/ECOR2_hypoxia_nfkb/"
library(magrittr)
ecor2.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hk.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-HK_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
hypoxia.pbs <- readr::read_csv(file = file.path(data.dir,"hypoxia_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)
ecor2i.pbs <- readr::read_csv(file = file.path(data.dir,"ECOR2-NFKBi_over_PBS.csv")) %>% dplyr::rename(SYMBOL = X1)

## add comparison ID
ecor2.pbs$comparison <- "E. coli"
hk.pbs$comparison <- "HK-E. coli"
hypoxia.pbs$comparison <- "hypoxia"
ecor2i.pbs$comparison <- "E. coli + SC-514"

## bind in single dataframe
data <- rbind(ecor2.pbs,
              hk.pbs,
              hypoxia.pbs,
              ecor2i.pbs)

plot.data <- data

## list of genes to extract
paths <- readr::read_csv("../results/ECOR2_hypoxia_nfkb/plotted-nfkb_complete-goANDreactome-results.csv")

adhesion <- paths[grep("Adherens junctions interactions", paths$Description),]$geneID

genes <- sapply(adhesion, 
                    function(x) {
                        v <- stringr::str_split(x, "/")
                                         })
genes <- c(unique(unlist(genes)))
genes <- genes[-grep("CDH9", genes)]

plot.data <- data[which(data$SYMBOL %in% genes),]

plot.data$category <- "Adherens junction"

adhesion <- paths[grep("cell-cell junction assembly", paths$Description),]$geneID
#adhesion <- paths[grep("cytoskeleton organization", paths$Description),]$geneID

genes <- sapply(adhesion, 
                    function(x) {
                        v <- stringr::str_split(x, "/")
                                         })
genes <- c(unique(unlist(genes)))
#genes <- genes[-grep("CDH9", genes)]

plot.data2 <- data[which(data$SYMBOL %in% genes),]

plot.data2$category <- "Cell-cell junction assembly"

plot.data <- rbind(plot.data, plot.data2)

library(magrittr)

## claculate zscore
scale_this <- function(x) as.vector(scale(x))

data.scaled <- dplyr::group_by(plot.data, SYMBOL) %>%
    dplyr::mutate(zscore = scale_this(log2FoldChange)) %>%
    dplyr::group_by(SYMBOL, comparison, category) %>%
    dplyr::summarize(mean_zscore = mean(zscore))

#data.scaled$category <- "Adherens junction interactions"
#data.scaled[grep("MUC", data.scaled$SYMBOL),]$category <- "Mucins"

data.scaled$comparison <- factor(data.scaled$comparison,
                                 levels = c(
                                            "E. coli", 
                                            "HK-E. coli",
                                            "E. coli + SC-514",
                                            "hypoxia"))

data.scaled <- data.scaled[order(data.scaled$comparison,data.scaled$mean_zscore),]
data.scaled$SYMBOL <- factor(data.scaled$SYMBOL,
                         levels = unique(data.scaled$SYMBOL))

## plot
library(ggplot2)
source("ggplot2-themes.R")

figure8a <- ggplot(data.scaled,
              aes(y = SYMBOL, x = comparison)) +
    geom_tile(stat = "identity", aes(fill = mean_zscore)) +
    facet_grid(category ~ ., scales = "free_y", space = "free", switch = "y", labeller = labeller(groupwrap = label_wrap_gen(10))) +
    scale_fill_distiller(name = "Z-score ", palette = "RdYlBu", breaks = c(-1, 0, 1)) +
    scale_y_discrete(position = "right") +
    ylab("") + xlab("") + 
    theme1 + 
    theme(strip.text =  element_text(size = 30),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 20),
          legend.position = "bottom",
	  legend.title = element_text(size = 18),
	  legend.key.size = unit(1,"cm"),
	  panel.spacing = unit(2, "lines"),
	  panel.border = element_blank()) +
    ggtitle("A")

png(filename = "../figures/figure8/figure8a.png", width = 400, height = 1200)
print(figure8a)
dev.off()
ggsave(filename = "../figures/figure8/eps/figure8a.eps", 
       plot = figure8a, 
       width = 8, height = 24)

## FIGURE 8 --------------------------------------------------------------------
## Figure 8B: FD4 permeability is reduced in /E. coli/ colonized HIOs treated with TNF and IFN 
## import data
## read in the raw data file from ImageJ
data <- readr::read_csv(file = "../data/figure8/results.csv")

## data wrangling ##########

## extract sample number from Label
sno <- "[_]([0-9]{2})[_R3D]"
## extract time point (image number) from label
time <- "t([0-9]{3})"
data$organoid <- as.numeric(stringr::str_match(data$Label, sno)[,2])
data$image <- as.numeric(stringr::str_match(data$Label, time)[,2])
## calculate minutes and hr timepoints
data$min <- (data$image * 10) - 10
data$hr <- data$min/60

## convert NaN to 0
is.nan.data.frame <- function(x) {
    do.call(cbind, lapply(x, is.nan))}
data$Mean[is.nan.data.frame(data$Mean)] <- 0

data$Mean_area <- data$Mean * data$Area

## subset baseline measurement
baseline <- subset(data, data$image == 1)
baseline$t0 <- baseline$Mean
baseline$t0_area <- baseline$Mean_area

baseline <- dplyr::select(baseline, organoid, t0, t0_area)
## create column with baseline measurements
data <- dplyr::left_join(data, baseline, by = "organoid") 

## read in sample key
key <- readr::read_csv(file = "../data/figure8/sample_key.csv")
## ... and merge
data <- dplyr::left_join(data, key, by = "organoid") 

## normalized fluorescence
data$normalized <- data$Mean/data$t0
data$normalized_area <- data$Mean_area/data$t0_area

library(magrittr)
data_mean <- dplyr::group_by(data, treatment, hr) %>%
    dplyr::summarise(mean = mean(normalized), 
                     stdev = sd(normalized),
                     num = n(),
                     iqr = IQR(normalized),
                     min = min(normalized),
                     max = max(normalized),
                     median = median(normalized))

data_mean$sem <- data_mean$stdev/sqrt(data_mean$num)


##Import csv
data <- read.table("../data/figure8/ECOR2_TNF_IFN_permeability.csv",header=TRUE,sep=",", stringsAsFactors = FALSE)

library(ggplot2)

#list of sample columns
samples <- grep("HIO",colnames(data))

#calculate row mean
data$mean <- rowMeans(data[,samples],na.rm=TRUE)
#calculate row SEM (counts number of non NA samples)
data$sem <- apply(data[,samples],1,sd,na.rm=TRUE)/sqrt(rowSums(!is.na(data[,samples])))
data$n <- rowSums(!is.na(data[,samples]))

## get the two datasets to play nicely
data <- dplyr::select(data, hr, mean, sem, treatment)
data_mean <- dplyr::select(data_mean, hr, mean, sem, treatment)
attributes(data_mean) <- NULL
data_mean <- rbind(data, data_mean)
data_mean <- data_mean[order(data_mean$treatment,data_mean$mean),]
## read in the raw data file from ImageJ
data <- readr::read_csv(file = "../data/figure8/041417_deltavision/thresh40.csv")

## data wrangling ##########

## extract sample number from Label
sno <- "[_]([0-9]{2})[_R3D]"
## extract time point (image number) from label
time <- "t([0-9]{3})"
data$organoid <- as.numeric(stringr::str_match(data$Label, sno)[,2])
data$image <- as.numeric(stringr::str_match(data$Label, time)[,2])
## calculate minutes and hr timepoints
data$min <- (data$image * 10) - 10
data$hr <- data$min/60

## convert NaN to 0
is.nan.data.frame <- function(x) {
    do.call(cbind, lapply(x, is.nan))}
data$Mean[is.nan.data.frame(data$Mean)] <- 0

data$Mean_area <- data$Mean * data$Area

## subset baseline measurement
baseline <- subset(data, data$image == 1)
baseline$t0 <- baseline$Mean
baseline$t0_area <- baseline$Mean_area

baseline <- dplyr::select(baseline, organoid, t0, t0_area)
## create column with baseline measurements
data <- dplyr::left_join(data, baseline, by = "organoid") 

## read in sample key
key <- readr::read_csv(file = "../data/figure8/041417_deltavision/sample_key.csv")
## ... and merge
data <- dplyr::left_join(data, key, by = "organoid") 

## normalized fluorescence
data$normalized <- data$Mean/data$t0
data$normalized_area <- data$Mean_area/data$t0_area

#data <- data[data$exclude != 1,]
#data <- subset(data, data$organoid > 3 & data$organoid < 15)
#data <- subset(data, data$organoid < 11 & data$organoid > 3)
library(magrittr)
data_mean2 <- dplyr::group_by(data, treatment, hr) %>%
    dplyr::summarise(mean = mean(normalized, na.rm = TRUE), 
                     stdev = sd(normalized, na.rm = TRUE),
                     num = n(),
                     iqr = IQR(normalized, na.rm = TRUE),
                     min = min(normalized, na.rm = TRUE),
                     max = max(normalized, na.rm = TRUE),
                     median = median(normalized, na.rm = TRUE))

data_mean2$sem <- data_mean2$stdev/sqrt(data_mean2$num)

data_mean2 <- dplyr::select(data_mean2, hr, mean, sem, treatment)
attributes(data_mean2) <- NULL
data_mean <- rbind(data_mean, data_mean2)

## read in the raw data file from ImageJ
data <- readr::read_csv(file = "../data/figure8/041717_deltavision/thesh40.csv")

## data wrangling ##########

## extract sample number from Label
sno <- "[_]([0-9]{2})[_R3D]"
## extract time point (image number) from label
time <- "t([0-9]{3})"
data$organoid <- as.numeric(stringr::str_match(data$Label, sno)[,2])
data$image <- as.numeric(stringr::str_match(data$Label, time)[,2])
## calculate minutes and hr timepoints
data$min <- (data$image * 10) - 10
data$hr <- data$min/60

## convert NaN to 0
is.nan.data.frame <- function(x) {
    do.call(cbind, lapply(x, is.nan))}
data$Mean[is.nan.data.frame(data$Mean)] <- 0

data$Mean_area <- data$Mean * data$Area

## subset baseline measurement
baseline <- subset(data, data$image == 1)
baseline$t0 <- baseline$Mean
baseline$t0_area <- baseline$Mean_area

baseline <- dplyr::select(baseline, organoid, t0, t0_area)
## create column with baseline measurements
data <- dplyr::left_join(data, baseline, by = "organoid") 

## read in sample key
key <- readr::read_csv(file = "../data/figure8/041717_deltavision/sample_key.csv")
## ... and merge
data <- dplyr::left_join(data, key, by = "organoid") 

## normalized fluorescence
data$normalized <- data$Mean/data$t0
data$normalized_area <- data$Mean_area/data$t0_area

## remove unreliable data
data <- subset(data, data$treatment != "PBS" & data$organoid != 7) ## HIOs dried out

library(magrittr)
data_mean3 <- dplyr::group_by(data, treatment, hr) %>%
    dplyr::summarise(mean = mean(normalized, na.rm = TRUE), 
                     stdev = sd(normalized, na.rm = TRUE),
                     num = n(),
                     iqr = IQR(normalized, na.rm = TRUE),
                     min = min(normalized, na.rm = TRUE),
                     max = max(normalized, na.rm = TRUE),
                     median = median(normalized, na.rm = TRUE))

data_mean3$sem <- data_mean3$stdev/sqrt(data_mean3$num)

data_mean3 <- dplyr::select(data_mean3, hr, mean, sem, treatment)
attributes(data_mean3) <- NULL
data_mean <- rbind(data_mean, data_mean3)

library(ggplot2)
source("ggplot2-themes.R")
data_mean$treatment <- factor(data_mean$treatment,
                              levels = c(
                                  "PBS",
				  "PBS + SC-514",
                                  "E. coli",
                                  "E. coli + TNF&IFN",
                                  "TNF&IFN",
                                  "E. coli + SC-514",
                                  "E. coli + TNF&IFN + SC-514",
                                  "LPS",
                                  "LPS + SC-514"))

plot.data1 <- subset(data_mean,
                     data_mean$treatment %in% c("PBS",
                                                "E. coli",
                                                "TNF&IFN",
                                                "E. coli + TNF&IFN"))

plot.data2 <- subset(data_mean,
                     data_mean$treatment %in% c("PBS",
                                                "E. coli",
                                                "PBS + SC-514",
                                                "E. coli + SC-514"))


figure8d <- ggplot(plot.data1, aes(x = hr, y = mean, color = treatment)) +
    scale_color_brewer(palette = "Spectral", name = "Treatment group", direction = -1) +
    geom_smooth(size = 3) +
    ylim(c(0.25, 1.05)) +
    xlab("Time post-treatment (h)") +
    ylab("Relative Fluorescence") +
    scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24)) +
    theme1 + ggtitle("D") +
    theme(
       # legend.position = "bottom",
          legend.position = c(0.3,0.13),
    	  legend.key.size = unit(1,"cm"),
          legend.text = element_text(size = 32))

png(filename = "../figures/figure8/figure8d.png", width = 800, height = 800)
print(figure8d)
dev.off()
ggsave(filename = "../figures/figure8/eps/figure8d.eps", 
       plot = figure8d, 
       width = 16, height = 16)

library(ggplot2)
source("ggplot2-themes.R")
data_mean$treatment <- factor(data_mean$treatment,
                              levels = c(
                                  "PBS",
                                  "E. coli",
                                  "E. coli + TNF&IFN",
                                  "PBS + SC-514",
                                  "TNF&IFN",
                                  "E. coli + SC-514",
                                  "E. coli + TNF&IFN + SC-514",
                                  "LPS",
                                  "LPS + SC-514"))

plot.data1 <- subset(data_mean,
                     data_mean$treatment %in% c("PBS",
                                                "E. coli",
                                                "TNF&IFN",
                                                "E. coli + TNF&IFN"))

plot.data2 <- subset(data_mean,
                     data_mean$treatment %in% c("PBS",
                                                "E. coli",
                                                "PBS + SC-514",
                                                "E. coli + SC-514"))



figure8b <- ggplot(plot.data2, aes(x = hr, y = mean, color = treatment)) +
    scale_color_brewer(palette = "Spectral", name = "Treatment group", direction = -1) +
    geom_smooth(size = 3) +
    ylim(c(0.25, 1.05)) +
    xlab("Time post-treatment (h)") +
    ylab("Relative Fluorescence") +
    scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24)) +
    theme1 + ggtitle("B") +
    theme(
       # legend.position = "bottom",
          legend.position = c(0.25,0.15),
    	  legend.key.size = unit(1,"cm"),
          legend.text = element_text(size = 32))

png(filename = "../figures/figure8/figure8b.png", width = 800, height = 800)
print(figure8b)
dev.off()
ggsave(filename = "../figures/figure8/eps/figure8b.eps", 
       plot = figure8b, 
       width = 16, height = 16)

## FIGURE 8 --------------------------------------------------------------------
## Figure 8C: Survival is impaired in presence of NFkB inhibitor
## import data
data <- readr::read_csv(file = "../data/figure1/161206_survival/survival_and_ELISA.csv")
## subset data
#data <- subset(data, data$plate > 3)
## create unique IDs
data$ID <- paste(data$well,data$plate, sep = "P")

## reshape dataframe
data <- dplyr::select(data, ID, day, bd1, bd2, il6, il8, vegf, treatment, Survival = dead)
data <- dplyr::rename(data, BD1 = bd1, BD2 = bd2, IL6 = il6, IL8 = il8, VEGF = vegf)
data <- reshape2::melt(data, id.vars = c('ID', 'day', 'treatment'))

## convert from pg/ml to pg/HIO
#data$value <- as.numeric(data$value, na.rm = TRUE)/2 # pg/1ml = 2*pg/0.5ml

## dotplot +/- SEM
library(magrittr)
data2 <- dplyr::group_by(data,day, treatment, variable) %>%
    dplyr::summarise(avg = mean(value,na.rm = TRUE),
              sem = sd(value, na.rm = TRUE)/n(), 
              total = sum(value, na.rm = TRUE))

data2[data2$variable == "Survival",]$avg <- 1 - data2[data2$variable == "Survival",]$avg
data2[data2$variable == "Survival",]$total <- 48 - data2[data2$variable == "Survival",]$total

levels(data2$variable) <- c("BD1 (pg/ml)", "BD2 (pg/ml)", "IL-6 (pg/ml)", "IL-8 (pg/ml)", "VEGF (pg/ml)", "Survival")

## list PBS first
data2 <- data2[data2$treatment !="HK- E. coli" ,]
data2 <- data2[complete.cases(data2),]
data2$treatment <- factor(data2$treatment,
                          levels = c("PBS", "E. coli", "E. coli + vehicle", "E. coli + SC-514"))

## survival numbers to reference in text
d9surv <- data2[data2$variable == "Survival" & data2$treatment == "E. coli" & data2$day == 9,]$total
d0surv <- data2[data2$variable == "Survival" & data2$treatment == "E. coli" & data2$day == 0,]$total
d3surv <- data2[data2$variable == "Survival" & data2$treatment == "E. coli" & data2$day == 3,]$total
d2surv <- data2[data2$variable == "Survival" & data2$treatment == "E. coli" & data2$day == 2,]$total

library(ggplot2)
source("ggplot2-themes.R")

figure8c  <- ggplot(data = data2[data2$variable == "Survival" & data2$treatment != "HK- E. coli",], 
                    aes(x = day, y = 1-avg, fill = treatment)) +
    geom_step(data = data2[data2$variable == "Survival" & data2$treatment != "HK- E. coli",],
              direction = "hv", aes(color = treatment), 
              size = 5) +
    xlab("Days post-microinjection") +
    ylab("Incidence of bacterial translocation") +
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7)) +
    ylim(c(0,1)) +
    scale_colour_brewer(palette = "Spectral", direction = -1) +
    guides(fill = guide_legend(title = NULL),
           color = guide_legend(title = NULL)) +
    theme1 + 
    #coord_fixed(ratio = 4) + 
    ggtitle("C") +
    theme(legend.position = c(0.2,0.85),
          legend.key.size = unit(1,"cm"),
	  legend.text = element_text(size = 32))

png(filename = "../figures/figure8/figure8c.png", width = 1000, height = 800)
print(figure8c)
dev.off()

ggsave(filename = "../figures/figure8/eps/figure8c.eps", 
       plot = figure8c, 
       width = 20, height = 16)

## Figure 7 multipanel ---------------------------------------------------------
source("ggplot2-themes.R")
source("custom_fun.R")
library(ggplot2)
library(gridExtra)

figure8e <- png2ggplot("../figures/figure8/figure8e.png") +
    img.theme + coord_fixed(ratio = 1/1.6) + ggtitle("E")

layout <- rbind(c(1,2,2,3,3,3),
                c(1,4,4,5,5,5))

## PDF output
pdf(file = "../figures/figure8/figure8_multipanel.pdf", width = 8250/300, height = 5250/300, onefile = FALSE)
gridExtra::grid.arrange(figure8a, figure8b, figure8c, figure8d, figure8e, layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/figure8/eps/figure8_multipanel.eps", 
       plot = gridExtra::grid.arrange(figure8a, figure8b, figure8c, figure8d, figure8e, layout_matrix = layout), 
       width = 44, height = 28)

## PNG output
png(filename = "../figures/figure8/figure8_multipanel.png", width = 2200, height = 1400)
gridExtra::grid.arrange(figure8a, figure8b, figure8c, figure8d,figure8e, layout_matrix = layout)
dev.off()

## FIGURE 8 --------------------------------------------------------------------
## Figure 8B: FD4 permeability is reduced in /E. coli/ colonized HIOs treated with TNF and IFN 
## import data
## read in the raw data file from ImageJ
data <- readr::read_csv(file = "../data/figure8/results.csv")

## data wrangling ##########

## extract sample number from Label
sno <- "[_]([0-9]{2})[_R3D]"
## extract time point (image number) from label
time <- "t([0-9]{3})"
data$organoid <- as.numeric(stringr::str_match(data$Label, sno)[,2])
data$image <- as.numeric(stringr::str_match(data$Label, time)[,2])
## calculate minutes and hr timepoints
data$min <- (data$image * 10) - 10
data$hr <- data$min/60

## convert NaN to 0
is.nan.data.frame <- function(x) {
    do.call(cbind, lapply(x, is.nan))}
data$Mean[is.nan.data.frame(data$Mean)] <- 0

data$Mean_area <- data$Mean * data$Area

## subset baseline measurement
baseline <- subset(data, data$image == 1)
baseline$t0 <- baseline$Mean
baseline$t0_area <- baseline$Mean_area

baseline <- dplyr::select(baseline, organoid, t0, t0_area)
## create column with baseline measurements
data <- dplyr::left_join(data, baseline, by = "organoid") 

## read in sample key
key <- readr::read_csv(file = "../data/figure8/sample_key.csv")
## ... and merge
data <- dplyr::left_join(data, key, by = "organoid") 

## normalized fluorescence
data$normalized <- data$Mean/data$t0
data$normalized_area <- data$Mean_area/data$t0_area

library(magrittr)
data_mean <- dplyr::group_by(data, treatment, hr) %>%
    dplyr::summarise(mean = mean(normalized), 
                     stdev = sd(normalized),
                     num = n(),
                     iqr = IQR(normalized),
                     min = min(normalized),
                     max = max(normalized),
                     median = median(normalized))

data_mean$sem <- data_mean$stdev/sqrt(data_mean$num)

data_mean_area <- dplyr::group_by(data, treatment, hr) %>%
    dplyr::summarise(mean = mean(normalized_area), 
                     stdev = sd(normalized_area),
                     num = n(),
                     iqr = IQR(normalized_area),
                     min = min(normalized_area),
                     max = max(normalized_area),
                     median = median(normalized_area))

data_mean_area$sem <- data_mean_area$stdev/sqrt(data_mean_area$num)

##Import csv
data <- read.table("../data/figure8/ECOR2_TNF_IFN_permeability.csv",header=TRUE,sep=",")

library(ggplot2)

#list of sample columns
samples <- grep("HIO",colnames(data))

#calculate row mean
data$mean <- rowMeans(data[,samples],na.rm=TRUE)
#calculate row SEM (counts number of non NA samples)
data$sem <- apply(data[,samples],1,sd,na.rm=TRUE)/sqrt(rowSums(!is.na(data[,samples])))
data$n <- rowSums(!is.na(data[,samples]))

## get the two datasets to play nicely
data$treatment <- gsub("E. coli", "ECOR2", data$treatment)
data <- dplyr::select(data, hr, mean, sem, treatment)
data_mean <- dplyr::select(data_mean, hr, mean, sem, treatment)
attributes(data_mean) <- NULL
data_mean <- rbind(data, data_mean)
data_mean <- data_mean[order(data_mean$treatment,data_mean$mean),]

library(ggplot2)
source("ggplot2-themes.R")
figure8b <- ggplot(data_mean, aes(x = hr, y = mean, color = treatment)) +
    scale_color_brewer(palette = "Spectral", name = "Treatment group", direction = -1) +
    geom_smooth(size = 3) +
    theme_bw() +
    ylim(c(0.3, 1.05)) +
    xlab("Time post-treatment (h)") +
    ylab("Relative Fluorescence") +
    scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24)) +
    theme1 + ggtitle("B") +
    theme(legend.position = c(0.25,0.2),
    	  legend.key.size = unit(1.5,"cm"),
          legend.text = element_text(size = 32))

png(filename = "../figures/figure8/figure8b.png", width = 800, height = 800)
print(figure8b)
dev.off()
ggsave(filename = "../figures/figure8/eps/figure8b.eps", 
       plot = figure8b, 
       width = 16, height = 16)
