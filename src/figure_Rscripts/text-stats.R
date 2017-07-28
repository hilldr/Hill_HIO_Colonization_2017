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
source("custom_fun.R")
stats <- lm_eqn(data[data$inject > 0 & data$fold !=0,],
                data[data$inject > 0 & data$fold !=0,]$inject,
                data[data$inject > 0 & data$fold !=0,]$fold)
  
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

## import data
data <- read.csv(file = "../data/figure1/ECOR2growth_fig1_timecourse.csv",
                 header =  TRUE, stringsAsFactors = FALSE)

## calculate means
group <- aggregate(CFU ~ hr, data = data, FUN = mean)
group.sem <- aggregate(CFU ~ hr, data = data, FUN = function(x) sd(x)/sqrt(length(x)))
group$sem <- group.sem$CFU

## stats
fit2d <- aov(CFU ~ hr, data = data[data$hr > 0,])

test2d.1 <- t.test(data[data$hr == 24 & data$CFU > 0,]$CFU,
                   data[data$hr == 0 & data$CFU > 0,]$CFU)

test2d.2 <- t.test(data[data$hr == 72 & data$CFU > 0,]$CFU,
                   data[data$hr == 24 & data$CFU > 0,]$CFU)
## Figure Figure 1F: /E. coli/ colonized HIOs are stable for up to 9 days 
## import data
data <- readr::read_csv(file = "../data/figure1/161206_survival/survival_and_ELISA.csv")
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
data2$treatment <- factor(data2$treatment,
                          levels = c("PBS", "E. coli"))

## survival numbers to reference in text
d9surv <- data2[data2$variable == "Survival" & data2$treatment == "E. coli" & data2$day == 9,]$total
d0surv <- data2[data2$variable == "Survival" & data2$treatment == "E. coli" & data2$day == 0,]$total
d3surv <- data2[data2$variable == "Survival" & data2$treatment == "E. coli" & data2$day == 3,]$total
d2surv <- data2[data2$variable == "Survival" & data2$treatment == "E. coli" & data2$day == 2,]$total

library(RColorBrewer)
library(ggplot2)
data <- read.table("../data/figure3/HIO_O2_final_data.csv",
                   header = TRUE, sep = ",", stringsAsFactors=FALSE)

## Summary stats
data.mean.0 <- aggregate(o2 ~ group, data[data$time == 0,], FUN = mean)
data.mean.0$sem <- aggregate(o2 ~ group, data[data$time == 0,],
                             FUN = function(x) sd(x)/sqrt(length(x)))$o2
data.mean.0$time <- 0

data.mean.24 <- aggregate(o2 ~ group, data[data$time == 24,],
                          FUN = mean)
data.mean.24$sem <- aggregate(o2 ~ group, data[data$time == 24,],
                              FUN = function(x) sd(x)/sqrt(length(x)))$o2
data.mean.24$time <- 24

data.mean.48 <- aggregate(o2 ~ group, data[data$time == 48,], FUN = mean)
data.mean.48$sem <- aggregate(o2 ~ group, data[data$time == 48,],
                              FUN = function(x) sd(x)/sqrt(length(x)))$o2
data.mean.48$time <- 48

t0 <- t.test(data[data$time == 0 & data$group == "PBS",]$o2,
             data[data$time == 0 & data$group != "PBS",]$o2,
             alternative = "greater")$p.value
t24 <- t.test(data[data$time == 24 & data$group == "PBS",]$o2,
              data[data$time == 24 & data$group != "PBS",]$o2,
              alternative = "greater")$p.value
t48 <- t.test(data[data$time == 48 & data$group == "PBS",]$o2,
              data[data$time == 48 & data$group != "PBS",]$o2,
              alternative = "greater")$p.value

## O2 concentrations in external media
media <- c(17.9, 18.4, 20.1, 19.2, 18.7)
media.mean <- mean(media)
media.sem <- sd(media)/sqrt(length(media))

base <- t.test(data[data$time == 0,]$o2,
             media,
             alternative = "less")$p.value

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

tt1.1 <- t.test2(m1 = data.sub[data.sub$comparison == 'E. coli',]$log2FoldChange[1],
               m2 = data.sub[data.sub$comparison == 'E. coli + SC-514',]$log2FoldChange[1],
               s1 = data.sub[data.sub$comparison == 'E. coli',]$lfcSE[1],
               s2 = data.sub[data.sub$comparison == 'E. coli + SC-514',]$lfcSE[1],
               n1 = 4, n2 = 4) 

tt2 <- t.test2(m1 = data.sub[data.sub$comparison == 'E. coli',]$log2FoldChange[1],
               m2 = data.sub[data.sub$comparison == 'hypoxia',]$log2FoldChange[1],
               s1 = data.sub[data.sub$comparison == 'E. coli',]$lfcSE[1],
               s2 = data.sub[data.sub$comparison == 'hypoxia',]$lfcSE[1],
               n1 = 4, n2 = 4)

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
hbd2.stats <- t.test(dat[dat$treatment == "E. coli",]$value,dat[dat$treatment == "PBS",]$value, alternative = "two.sided")
hbd2.stats2 <- t.test(dat[dat$treatment == "E. coli",]$value,dat[dat$treatment == "E. coli + SC-514",]$value, alternative = "two.sided")
hbd2.stats3 <- t.test(dat[dat$treatment == "E. coli",]$value,dat[dat$treatment == "heat-inactivated",]$value, alternative = "two.sided")
hbd2.stats4 <- t.test(dat[dat$treatment == "E. coli",]$value,dat[dat$treatment == "hypoxia",]$value, alternative = "two.sided")

#+begin_src R :session *R* :results silent :exports none :eval yes :tangle figure_Rscripts/figure5.R
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

library(magrittr) 
plot.data <- data2[data2$dose == 1e-8| data2$dose == 0.1 |data2$dose == 1,] %>%
    dplyr::group_by(dose, Time) %>%
    dplyr::summarize(avg = mean(value), sem = sd(value)/sqrt(n()))

plot.data[plot.data$dose == 1e-8,]$dose <- 0

#+begin_src R :session *R* :results silent :exports none :eval yes :tangle figure_Rscripts/figure7.R
## FIGURE 7 --------------------------------------------------------------------
## Figure 7B&D: FD4 permeability is reduced in /E. coli/ colonized HIOs treated with TNF and IFN 
## import data
## read in the raw data file from ImageJ
data <- readr::read_csv(file = "../data/figure7/results.csv")

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
key <- readr::read_csv(file = "../data/figure7/sample_key.csv")
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
data <- read.table("../data/figure7/ECOR2_TNF_IFN_permeability.csv",header=TRUE,sep=",", stringsAsFactors = FALSE)

library(ggplot2)

#list of sample columns
samples <- grep("HIO",colnames(data))

#calculate row mean
data$mean <- rowMeans(data[,samples],na.rm=TRUE)
#calculate row SEM (counts number of non NA samples)
data$stdev <- apply(data[,samples],1,sd,na.rm=TRUE)
data$sem <- data$stdev/sqrt(rowSums(!is.na(data[,samples])))
data$num <- rowSums(!is.na(data[,samples]))

## get the two datasets to play nicely
data <- dplyr::select(data, hr, mean, sem, stdev, num, treatment)
data_mean <- dplyr::select(data_mean, hr, mean, sem, stdev, num, treatment)
attributes(data_mean) <- NULL
data_mean <- rbind(data, data_mean)
data_mean <- data_mean[order(data_mean$treatment,data_mean$mean),]
## read in the raw data file from ImageJ
data <- readr::read_csv(file = "../data/figure7/041417_deltavision/thresh40.csv")

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
key <- readr::read_csv(file = "../data/figure7/041417_deltavision/sample_key.csv")
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

data_mean2 <- dplyr::select(data_mean2, hr, mean, sem, stdev, num, treatment)
attributes(data_mean2) <- NULL
data_mean <- rbind(data_mean, data_mean2)

## read in the raw data file from ImageJ
data <- readr::read_csv(file = "../data/figure7/041717_deltavision/thesh40.csv")

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
key <- readr::read_csv(file = "../data/figure7/041717_deltavision/sample_key.csv")
## ... and merge
data <- dplyr::left_join(data, key, by = "organoid") 

## normalized fluorescence
data$normalized <- data$Mean/data$t0
data$normalized_area <- data$Mean_area/data$t0_area

#data <- data[data$exclude != 1,]
#data <- subset(data, data$organoid > 3 & data$organoid < 15)
data <- subset(data, data$treatment != "PBS") ## HIOs dried out
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

data_mean3 <- dplyr::select(data_mean3, hr, mean, sem, stdev, num, treatment)
attributes(data_mean3) <- NULL
data_mean <- rbind(data_mean, data_mean3)

#T-tests
source("custom_fun.R")
## t-test
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
# t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
tt.fig7b <- t.test2(m1 = data_mean[data_mean$treatment == "E. coli" &
                                       data_mean$hr == 20,]$mean[1],
                    m2 = data_mean[data_mean$treatment == "E. coli + SC-514" &
                                       data_mean$hr == 20,]$mean,
                    s1 = data_mean[data_mean$treatment == "E. coli" &
                                       data_mean$hr == 20,]$stdev[1],
                    s2 = data_mean[data_mean$treatment == "E. coli + SC-514" &
                                       data_mean$hr == 20,]$stdev[1],
                    n1 = data_mean[data_mean$treatment == "E. coli" &
                                       data_mean$hr == 20,]$num[1],
                    n2 = data_mean[data_mean$treatment == "E. coli + SC-514" &
                                       data_mean$hr == 20,]$num,
                    equal.variance = TRUE)
tt.fig7d <- t.test2(m1 = data_mean[data_mean$treatment == "PBS" &
                                       data_mean$hr == 20,]$mean[2],
                    m2 = data_mean[data_mean$treatment == "TNF&IFN" &
                                       data_mean$hr == 20,]$mean[1],
                    s1 = data_mean[data_mean$treatment == "PBS" &
                                       data_mean$hr == 20,]$stdev[2],
                    s2 = data_mean[data_mean$treatment == "TNF&IFN" &
                                       data_mean$hr == 20,]$stdev[1],
                    n1 = data_mean[data_mean$treatment == "PBS" &
                                       data_mean$hr == 20,]$num[2],
                    n2 = data_mean[data_mean$treatment == "TNF&IFN" &
                                       data_mean$hr == 20,]$num[1],
                    equal.variance = TRUE)
tt.fig7d2 <- t.test2(m1 = data_mean[data_mean$treatment == "E. coli + TNF&IFN" &
                                       data_mean$hr == 20,]$mean[2],
                    m2 = data_mean[data_mean$treatment == "TNF&IFN" &
                                       data_mean$hr == 20,]$mean[1],
                    s1 = data_mean[data_mean$treatment == "E. coli + TNF&IFN" &
                                       data_mean$hr == 20,]$stdev[2],
                    s2 = data_mean[data_mean$treatment == "TNF&IFN" &
                                       data_mean$hr == 20,]$stdev[1],
                    n1 = data_mean[data_mean$treatment == "E. coli + TNF&IFN" &
                                       data_mean$hr == 20,]$num[2],
                    n2 = data_mean[data_mean$treatment == "TNF&IFN" &
                                       data_mean$hr == 20,]$num[1],
                    equal.variance = TRUE)

## FIGURE 7 --------------------------------------------------------------------
## Figure 7C: Survival is impaired in presence of NFkB inhibitor
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

data <- data[data$variable == "Survival",]
library(survival)
f7s <- Surv(time = data$day, event = data$value)
model <- coxph(f7s ~ treatment, data = data)

library(magrittr)
data2 <- dplyr::group_by(data,day, treatment, variable) %>%
    dplyr::summarise(avg = mean(value,na.rm = TRUE),
              sem = sd(value, na.rm = TRUE)/n(), 
              total = sum(value, na.rm = TRUE))
