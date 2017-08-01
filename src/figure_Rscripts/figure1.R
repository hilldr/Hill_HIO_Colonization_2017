## FIGURE 1 --------------------------------------------------------------------
## Figure 1A: /E. coli/ density increases rapidly following HIO microinjection
## import data
data <- readr::read_csv(file = "../data/figure1/010716_01_R3D.csv")

## plot
library(ggplot2)
library(scales)
source("ggplot2-themes.R")

figure1a <- ggplot(data = data[data$hr <= 18,], aes(x = hr, y = Median*Area)) +
    geom_point(shape = 21, fill = color.set[3], size = 8) +
    xlab("Time post-microinjection (h)") +
    ylab(expression("Mean fluorescent intensity" %*% "pixel area")) +   
    scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15, 18)) +       
    scale_y_log10(
        breaks = trans_breaks("log10", function(x) 10^x)(c(1, 1e7)),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l", size = 1,
                        short = unit(.75,"mm"),
                        mid = unit(3,"mm"),
                        long = unit(5,"mm"))+       
    theme1 + ggtitle("A")

source("custom_fun.R")
fig <- png2ggplot("../data/figure1/010716/010716_01_R3D_w594_t01.png") +
    img.theme + coord_fixed(ratio = 1/4)

library(gridExtra)
layout <- rbind(c(2,2),
                c(2,2),
                c(2,2),
                c(1,1))


png(filename = "../figures/figure1/figure1a.png", width = 800, height = 800)

grid.arrange(fig, figure1a,
             layout_matrix = layout)
dev.off()

ggsave(filename = "../figures/figure1/eps/figure1a.eps", 
       plot = gridExtra::grid.arrange(fig, figure1a,
             layout_matrix = layout), 
       width = 16, height = 16)
         

## FIGURE 1C --------------------------------------------------------------------
## Figure 1C: Figure 1C: Maximum /E. coli/ density is 10^{5} CFU/HIO
## import data
## import data
data <- read.csv(file = "../data/figure1/ECOR2growth_fig1.csv",
                 header =  TRUE, stringsAsFactors = FALSE)

## Index of CFU/HIO injected for each Sample condition (A-G)
sample.table <- read.table("../data/figure1/sample_table_fig1.csv",
                           header = TRUE, sep = ",", stringsAsFactors = FALSE)

## Generate index of rows in sample table that match the sample labels in data
id <- match(data$sample,sample.table$sample)
## create column in data of of CFU/HIO values in sample table in matching rows listed in id
data$inject <- sample.table[id,]$value
data$fold <- data$mean/data$inject
data$increase <- ifelse(data$mean > data$inject,"increase","decrease")
## generate stats string for plot
source("custom_fun.R")
stats <- lm_eqn(data[data$inject > 0 & data$fold !=0,],
                data[data$inject > 0 & data$fold !=0,]$inject,
                data[data$inject > 0 & data$fold !=0,]$fold)

library(ggplot2)
library(scales)
source("ggplot2-themes.R")

figure1c <- ggplot(data, aes(x = inject, y = fold)) +
    geom_smooth(data = data[data$inject > 0 & data$fold !=0,],
                aes(x=inject, y=fold), colour = "black",
                size = 2,
                method = "lm",
                formula = y ~ x,
                level = 0.95) +
    geom_hline(yintercept = 1, color = "black", size = 0.5, linetype = "dashed") +
    geom_point(size = 8,shape = 21, fill=color.set[2]) + 
    ylab(latex2exp::TeX("$\\textbf{$\\frac{CFU$\\cdot{}HIO_{$\\textit{t}=24}^{-1}}{CFU$\\cdot{}HIO_{$\\textit{t}=0}^{-1}}}$")) +
    ggtitle("C") + 
    xlab("CFU injected per HIO") + theme1 + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(limits = c(10, 1e5),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides="bl",size= 2) +
    scale_fill_brewer(palette = "Set1") + 
    scale_color_brewer(palette = "Set1") +
    ## size of stats label
    annotate("text",x = 1000, y = 100000,
                           label = substr(stats,62,150),
                           parse = TRUE,
                           size = 10) +
    theme(axis.title.y = element_text(vjust =-0.8))

png(filename = "../figures/figure1/figure1c.png", width = 800, height = 800)
print(figure1c)
dev.off()

ggsave(filename = "../figures/figure1/eps/figure1c.eps", 
       plot = figure1c, 
       width = 16, height = 16)

## FIGURE 1D --------------------------------------------------------------------
## Figure 1D: E. coli grows within the HIO lumen
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

## plot Figure 1D
library(ggplot2)
library(grid)
library(scales)
source("ggplot2-themes.R")

figure1d <- ggplot(data[data$CFU > 0,], aes(x = hr, y = CFU)) +
    geom_boxplot(size = 2, aes(group = hr), fill = color.set[2]) +
    geom_point(size = 8, shape = 21, fill = color.set[2]) +
    ylab(latex2exp::TeX("$\\textbf{CFU$\\cdot{}HIO^{-1}}$")) +
    ggtitle("D") + 
    xlab("Time post-microinjection (h)") + theme1 +
    scale_y_log10(limits = c(1,50000000),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(breaks = c(0, 24, 48, 72)) +
    annotation_logticks(sides = "l", size = 2,
                        short = unit(.75,"mm"),
                        mid = unit(3,"mm"),
                        long = unit(5,"mm")) +
    annotate(geom = "segment",
             x = 24, xend = 72, y = 10^7, yend = 10^7,
             lwd = 2.5) +
    annotate(geom = "text", x = 48, y = 10^7,
             label = paste("P =",format(test2d.2$p.value, digits = 2)), 
             size = 10, vjust = -0.5) + 
    scale_fill_brewer(palette = "Set1") + 
    scale_color_brewer(palette = "Set1")


png(filename = "../figures/figure1/figure1d.png", width = 800, height = 800)
print(figure1d)
dev.off()

ggsave(filename = "../figures/figure1/eps/figure1d.eps", 
       plot = figure1d, 
       width = 16, height = 16)

## FIGURE 1 --------------------------------------------------------------------
## Figure Figure 1F: /E. coli/ colonized HIOs are stable for up to 9 days 
## import data
data <- readr::read_csv(file = "../data/figure1/161206_survival/survival_and_ELISA.csv")
## subset data
data <- subset(data, data$treatment == "PBS" | data$treatment == "E. coli")
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

library(ggplot2)
source("ggplot2-themes.R")

figure1f  <- ggplot(data = data2[data2$variable == "Survival",], 
                    aes(x = day, y = 1-avg, fill = treatment)) +
    geom_step(data = data2[data2$variable == "Survival",],
              direction = "hv", aes(color = treatment), 
              size = 5) +
    xlab("Days post-microinjection") +
    ylab("Bacterial translocation rate") +
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9)) +
    ylim(c(0,1)) +
    scale_colour_brewer(palette = "Set1", direction = -1) +
    scale_fill_brewer(palette = "Set1") +
    guides(fill = guide_legend(title = NULL),
           color = guide_legend(title = NULL)) +
    theme1 + 
    #coord_fixed(ratio = 4) + 
    ggtitle("F") +
    theme(legend.position = c(0.2, 0.8),
          legend.key.size = unit(2,"cm"),
	  legend.text = element_text(size = 32))

png(filename = "../figures/figure1/figure1f.png", width = 1000, height = 800)
print(figure1f)
dev.off()

ggsave(filename = "../figures/figure1/eps/figure1f.eps", 
       plot = figure1f, 
       width = 20, height = 16)

## Figure 1 multipanel ---------------------------------------------------------
library(ggplot2)
library(gridExtra)
source("ggplot2-themes.R")
source("custom_fun.R")

figure1b <- png2ggplot("../figures/figure1/figure1b.png") +
    img.theme + ggtitle("B")

figure1e <- png2ggplot("../figures/figure1/figure1e.png") +
    img.theme + ggtitle("E") + coord_fixed(ratio = 0.6)

layout2 <- rbind(c(2,2),
                c(2,2),
                c(2,2),
                c(1,1))

layout <- rbind(c(1,1,1,2,2,2),
                c(1,1,1,2,2,2),
                c(3,3,4,4,5,5),
                c(3,3,4,4,6,6))

figure1a <- grid.arrange(fig, figure1a,
                         layout_matrix = layout2)

## PDF output
pdf(file = "../figures/figure1/figure1_multipanel.pdf", width = 7500/300, height = 7500/300, onefile = FALSE)
gridExtra::grid.arrange(figure1a, figure1b, figure1c, figure1d, figure1e, figure1f,
             layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/figure1/eps/figure1_multipanel.eps", 
       plot = gridExtra::grid.arrange(figure1a, figure1b, figure1c, figure1d, figure1e, figure1f, layout_matrix = layout), 
       width = 20, height = 20)

## PNG output
png(filename = "../figures/figure1/figure1_multipanel.png", width = 2000, height = 2000)
gridExtra::grid.arrange(figure1a, figure1b, figure1c, figure1d, figure1e, figure1f,
             layout_matrix = layout)
dev.off()
