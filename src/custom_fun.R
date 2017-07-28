## Revised GSEA plot function --------------------------------------------------
## modified from https://github.com/GuangchuangYu/clusterProfiler
gseaplot2 <- function (gseaResult, geneSetID, by = "all") 
{
    by <- match.arg(by, c("runningScore", "position", "all"))
    x <- ymin <- ymax <- runningScore <- es <- pos <- geneList <- NULL
    p <- ggplot(gseaResult, geneSetID = geneSetID,
                aes(x = x, 
                    ymin = ymin, ymax = ymax))  +
        xlab("Position in the Ranked List of Genes")
    if (by == "runningScore" || by == "all") {
        p.res <- p 
        p.res <- p.res + geom_area(aes(y = runningScore), fill = "grey90")
        p.res <- p.res + geom_line(aes(y = runningScore), size = 3)
	p.res <- p.res + geom_linerange(colour = color.set[1])
        enrichmentScore <- gseaResult@result[geneSetID, "enrichmentScore"]
        es.df <- data.frame(es = which(p$data$runningScore == 
            enrichmentScore))
        p.res <- p.res + geom_vline(data = es.df, aes(xintercept = es), 
            colour = color.set[1], linetype = "dashed")
        p.res <- p.res + ylab("Running Enrichment Score")
        p.res <- p.res + geom_hline(aes(yintercept = 0))
    }
    if (by == "position" || by == "all") {
        df2 <- data.frame(pos = which(p$data$position == 1))
        p.pos <- p + geom_vline(data = df2, aes(xintercept = pos), 
            colour = "#DAB546", alpha = 0.3)
        p.pos <- p.pos + geom_line(aes(y = geneList), colour = "red")
        p.pos <- p.pos + ylab("Phenotype")
        p.pos <- p.pos + geom_hline(aes(yintercept = 0))
    }
    if (by == "runningScore") 
        return(p.res)
    if (by == "position") 
        return(p.pos)

    p.res <- p.res + theme(axis.title.x = element_text(vjust = -0.3))
    return(p.res)
}

## Convert PNG to ggplot object ------------------------------------------------
png2ggplot <- function(filename) {
    library(ggplot2)
    img <- png::readPNG(filename, native = TRUE)
    grob <- grid::rasterGrob(img, interpolate = FALSE)
    fig <- qplot(1:100, 1:100, alpha = I(0)) +
        theme_bw() +
      #  geom_point(size = 0, color = "white") +
        annotation_custom(grob, xmin = -Inf,
                          xmax = Inf,
                          ymin = -Inf,
                          ymax = Inf) 
        #coord_fixed(ratio = 1) 
    return(fig)
}

## GET EQUATION AND R-SQUARED AS STRING ----------------------------------------
## SOURCE: http://goo.gl/K4yh
lm_eqn <- function(df,x,y){
    m <- lm(log(y) ~ log(x), df)
    eq <- substitute(italic(y) == 10^(a + b %.% italic(log(x)))*","~~italic(r)^2~"="~r2*","~~italic(P)~"="~pv, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2),
                          pv = format(anova(m)$'Pr(>F)'[1],digits = 3), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq))                 
}

lm_eqn2 <- function(df,x,y){
    m <- lm(y ~ log(x), df, na.action = 'na.exclude')
    eq <- substitute(italic(y) == a + b %.% italic(log(x))*","~~italic(r)^2~"="~r2*","~~italic(P)~"="~pv, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2),
                          pv = format(anova(m)$'Pr(>F)'[1],digits = 3), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq))                 
}

## function to format decimals as precentage -----------------------------------
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
## t-test
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE ) 
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
        df <- n1+n2-2
    }      
    t <- (m1-m2-m0)/se 
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat) 
}

## create placeholder ggplot with test
fillerggplot <- function(text) {
    library(ggplot2)
    fig <- qplot(1:100, 1:100, alpha = I(0)) +
        theme_bw() + img.theme +
        annotate("text", x = 50, y = 50, label = text, size = 20)
      #  geom_point(size = 0, color = "white") +
        #coord_fixed(ratio = 1) 
    return(fig)
}
