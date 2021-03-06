## FIGURE 3 --------------------------------------------------------------------
## Figure 4A: /E. coli/ colonization is associated with a reduction in luminal O_{2}
## import data
library(RColorBrewer)
library(ggplot2)
data <- read.table("../data/figure4/HIO_O2_final_data.csv",
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

summary.data <- rbind(data.mean.0,data.mean.24,data.mean.48)

t0 <- t.test(data[data$time == 0 & data$group == "PBS",]$o2,
             data[data$time == 0 & data$group != "PBS",]$o2,
             alternative = "greater")$p.value
t24 <- t.test(data[data$time == 24 & data$group == "PBS",]$o2,
              data[data$time == 24 & data$group != "PBS",]$o2,
              alternative = "greater")$p.value
t48 <- t.test(data[data$time == 48 & data$group == "PBS",]$o2,
              data[data$time == 48 & data$group != "PBS",]$o2,
              alternative = "greater")$p.value

stats <- cbind(summary.data[summary.data$group == "+ E. coli",],c(t0,t24,t48))
colnames(stats)[5] <- "p"
star <- function(x){
    if (x > 0.05){
        return("")
    } else {
        if (x < 0.001){
            return("***")
        } else {
            if (x < 0.01){
                return("**")
            } else {
                return("*")
            }
        }
    }
}

stats$star <- as.character(lapply(stats$p, star))

## plot
library(ggplot2)  
source("ggplot2-themes.R")

data$group <- factor(data$group,
                     levels = c("PBS", "+ E. coli"))


figure4a <- ggplot(data[data$time < 72,],
                   aes(x= factor(time), y = o2, fill = treatment)) +
    geom_boxplot(notch = FALSE, aes(fill = group),
                 outlier.colour = NULL, outlier.shape = 21, outlier.size = 5,
                 size = 2, width = 0.4) +
    ylab(expression("%O"[2])) +
    ggtitle("A") + 
    xlab("Hours post-microinjection") + 
    theme1 +
    theme(legend.position = c(0.18, 0.12),
          legend.text = element_text(size = 28)) +
    guides(fill = guide_legend(override.aes = list(size = 1))) +

    annotate("text", x = c(1,2,3), y = c(9.75,6,9.75),
             label = paste("P = ",round(stats$p, digits = c(2,3,6))), size = 10) +
    annotate("segment",
             x = 1.85, xend = 2.15, y = 5.5, yend = 5.5,
             color = "black", size = 2) +
    annotate("segment",
             x = 2.85, xend = 3.15, y = 9.25, yend = 9.25,
             color = "black", size = 2) + 
    scale_fill_brewer(palette = "Set1", direction = -1) +
    scale_color_brewer(palette = "Set1") +
    theme(legend.key.size = unit(2,"cm"),
	  legend.text = element_text(size = 32))

png(filename = "../figures/figure4/figure4a.png", width = 800, height = 800)
print(figure4a)
dev.off()

ggsave(filename = "../figures/figure4/eps/figure4a.eps", 
       plot = figure4a, 
       width = 20, height = 20)

## FIGURE 3 --------------------------------------------------------------------
## Figure 4B: Reduction in luminal O_{2} associated with /E. coli/ growth
## compute fit line for CFU ~ O2
source("custom_fun.R")
stats.b <- lm_eqn2(data[data$cfu !=0 & data$time != 0 & data$time != 72,],
                log(data[data$cfu !=0 & data$time != 0 & data$time != 72,]$cfu),
                data[data$cfu !=0 & data$time != 0 & data$time != 72,]$o2)

library(scales)
library(ggplot2)
source("ggplot2-themes.R")

figure4b <- ggplot(data[data$cfu != 0 & data$time != 0 & data$time != 72,],
                   aes(x = cfu, y = o2)) +
    geom_smooth(data = data[data$cfu !=0 & data$time != 0 & data$time != 72,],
                aes(x = cfu, y = o2), colour = "black", size = 2,
                method = "lm",
                formula = y ~ x,
                level = 0.95) +
    geom_point(size = 8,shape = 21, fill = color.set[1] ) +
    scale_y_continuous(breaks = c(0:5)) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "b", size = 2,
                        short = unit(.75,"mm"),
                        mid = unit(3,"mm"),
                        long = unit(5,"mm")) +
    ylab(expression("%O"[2])) +
    ggtitle("B") + 
    xlab("CFU/HIO") + 
    theme1 +
    annotate("text",x = 10000, y = 4.5, label = stats.b, parse = TRUE, size = 8)

png(filename = "../figures/figure4/figure4b.png", width = 1000, height = 1000)
print(figure4b)
dev.off()

ggsave(filename = "../figures/figure4/eps/figure4b.eps", 
       plot = figure4b, 
       width = 20, height = 20)

## PMDZ staining table
stain <- data.frame(pbs = c(1,12,13, NA),
                    ecoli = c(10,2,12, round(chisq.test(stain[1:2,c(1,2)])$p.value, digits = 3)),
                    hk = c(3,10,13, round(chisq.test(stain[1:2,c(1,3)])$p.value, digits = 3)),
                    o2 = c(14,0,14,round(chisq.test(stain[1:2,c(1,4)])$p.value, digits = 5)))

rownames(stain) <- c("PMDZ +", "PMDZ -", "N", "P")
colnames(stain) <- c("PBS (negative control)", "E. coli (live)", "E. coli (heat-inactivated)", "1% oxygen (positive control)")

library(gridExtra)
library(grid)
tt <- ttheme_default(
		     colhead=list(fg_params=list(col="black", fontface=4L, cex = 3),
                                  bg_params = list(fill = c(rep('grey70', times = 2), 'grey80', 'grey80'))),
                     core = list(fg_params = list(fontface = "bold", cex = 3)),
                     rowhead = list(fg_params=list(col="black", fontface=4L, cex = 3, hjust = 0.95),
                                    bg_params = list(fill = c("white",rep(c('grey90','grey80'), times = 5)))))
grid.table(t(stain), theme = tt) 

## Figure 4 multipanel ---------------------------------------------------------
source("ggplot2-themes.R")
library(ggplot2)
library(gridExtra)

figure4c <- png2ggplot("../figures/figure4/figure4c.png") +
    img.theme + ggtitle("C") + coord_fixed(ratio = 0.5)

layout <- rbind(c(1,1,2,2),
                c(1,1,2,2),
                c(3,3,3,3),
                c(3,3,3,3),
		c(4,4,4,4))

## PDF output
pdf(file = "../figures/figure4/figure4_multipanel.pdf", width = 6000/300, height = 7200/300, onefile = FALSE)
gridExtra::grid.arrange(figure4a, figure4b, figure4c, tableGrob(t(stain), theme = tt), 
layout_matrix = layout)
dev.off()

## EPS output
ggsave(filename = "../figures/figure4/eps/figure4_multipanel.eps", 
       plot = gridExtra::grid.arrange(figure4a, figure4b, figure4c, tableGrob(t(stain), theme = tt), layout_matrix = layout), 
       width = 20, height = 20)

## PNG output
png(filename = "../figures/figure4/figure4_multipanel.png", width = 1600, height = 1600)
gridExtra::grid.arrange(figure4a, figure4b, figure4c, tableGrob(t(stain), theme = tt), layout_matrix = layout)
dev.off()
