#' # SETUP

#' Load packages and paths
#+ setup-packages
library(grid)
library(ggplot2)
library(plyr)
library(RColorBrewer)


#' ## Color Scheme
#' 
#' CMI Recommendations
#+ setup-colors
cmi_main_blue="#0071b2"
cmi_grey="#929d9e"
cmi_light_blue="#00c4d9"
cmi_pea_green="#b5bf00"

cmi_rich_green="#73933d"
cmi_rich_purple="#8e7fac"
cmi_rich_red="#d75920"
cmi_rich_blue="#4c87a1"
cmi_rich_aqua="#66c7c3"
cmi_rich_orange="#eebf42"

cmi_vibrant_yellow="#ffd457"
cmi_vibrant_orange="#f58025"
cmi_vibrant_green="#78a22f"
cmi_vibrant_garnet="#e6006f"
cmi_vibrant_purple="#9A4d9e"
cmi_vibrant_blue="#19398a"

cmi_site_colors = c(cmi_vibrant_blue,
                    cmi_rich_blue,
                    cmi_vibrant_purple,
                    cmi_vibrant_garnet,
                    cmi_rich_red,
                    cmi_vibrant_orange,
                    cmi_vibrant_yellow,
                    cmi_vibrant_green)
cmi_site_colors_ramp = colorRampPalette(cmi_site_colors)


#' ## Load data
#'
#' Read in the data and then some
#+ setup-load
#setwd("/home2/zarrar/projects/qc")
#setwd("~/zarrar/qc")
setwd("~/Dropbox/Research/cmi/qc")
df      <- read.csv("corr.qc/qc_filt_epi_derivatives.csv")[,-1]
nsites  <- length(unique(df$site))

#' ## Percentiles
#'
#' In our plots, we want to have percentile lines for each QC measure to 
#' indicate the distribution of each site relative to the whole sample
#+ setup-percentile
qc.measures <- colnames(df)[!(colnames(df) %in% c("uniqueid", "subid", "site", "site.name", "session", "scan", "global"))]
qvals       <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
qcat        <- c(1,5,25,50,25,5,1)
qline       <- c(3, 2, 5, 1, 5, 2, 3)
qsize       <- c(.4, .25, .3, .25, .3, .25, .4)
qcols       <- c("grey10", "grey10", "grey10", "grey50", "grey10", "grey10", "grey10")
#qcols       <- brewer.pal(8, "Dark2")[c(1,2,3,4,3,2,1)]

#' Now let's get the percentiles
#+ percentile
percentiles <- apply(subset(df, select=qc.measures), 2, quantile, qvals, na.rm=TRUE)
percentiles <- as.data.frame(cbind(percentiles, qcat, qline, qsize))
percentiles$qline <- as.factor(qline)
percentiles$qcat  <- as.factor(qcat)
print(percentiles)

#' ## Measure Descriptions
#' Associate a detailed description with each measure
#+ setup-measure
dnames  <- c("fALFF", "REHO", "VMHC")
mnames  <- c("Median", "75th Percentile", "90th Percentile", "FWHM (mm)", "Mean")
descs   <- expand.grid(m=mnames, d=dnames) 
descs   <- paste(descs[,2], descs[,1], sep=" - ")
mdf <- data.frame(
    measure = qc.measures, 
    description = descs
)
print(mdf)
cat("CHECK ABOVE...do the columns for each row match\n")

#' # PLOTS

#' ## More Setup
#' 
#' This function will add percentile lines in the background
#' plot: ggplot object
#' pdf: percentile data frame
#+ plot-function
compile_percentiles <- function(pdf, measure, cols=NULL) {
  ret <- lapply(1:nrow(pdf), function(i) {
    p <- pdf[i,]
    if (!is.null(cols)) {
      plot <- geom_hline(aes_string(yintercept=measure), data=p, 
                         size=as.numeric(p$qsize), linetype=as.numeric(p$qline), color=cols[i])
      #as.character(p$qcolor[1])
    } else {
      plot <- geom_hline(aes_string(yintercept=measure), data=p, 
                         size=as.numeric(p$qsize[1]), linetype=as.numeric(p$qline[1]), color="grey50")
    }
    return(plot)
  })
  return(ret)
}


#' ## Outliers
#' Sometimes extreme data-points can skew the plot 
#' and make it difficult to see the spread of the data.
#'  
#' I will avoid plotting those outlier points by setting the axis
#' to only the points that I want. ggplot will complain but let her (poor baby).
#'
#+ plot-find-outliers
# functions
range.outlier.iqr <- function(x, times=3) {
    upper.limit <- quantile(x, 0.75) + times*IQR(x)
    lower.limit <- quantile(x, 0.25) - times*IQR(x)
    return(c(lower.limit, upper.limit))
}
outlier.iqr <- function(x, times=3) {
    tmp         <- range.outlier.iqr(x, times)
    lower.limit <- tmp[1]
    upper.limit <- tmp[2]
    return((x > upper.limit) | (x < lower.limit))
}
# outlier values (if any)
lst.outlier.iqr <- llply(qc.measures, function(measure) {
    ret     <- subset(df, select=c("uniqueid", "subid", "site", "site.name", "session", "scan", measure))
    inds    <- outlier.iqr(df[[measure]])
    return(ret[inds,])
})
names(lst.outlier.iqr) <- qc.measures
# new ranges of our plots (sans outliers)
df.range.iqr <- as.data.frame(sapply(qc.measures, function(m) {
    inds <- !outlier.iqr(df[[m]])
    range(df[[m]][inds]) * c(0.99,1.01)
}))


#' ## QC Derivative Measures
#' I will be plotting a bunch of different data-sets here. First I will have 
#' all the data, then I will have it remote 3 x the IQR. So two sets of the
#' same plots.
#'
#' ### REMOVE OUTLIERS
#+ plot-qc-no-outliers, fig.width=15, fig.height=10
for (i in 1:nrow(mdf)) {
  measure <- as.character(mdf$measure[i])
  desc <- as.character(mdf$description[i])
  
  # ### Option 1
  # First, I'll plot ones with 1%, 5%, 25%, and 50% percentile lines
  pg1=ggplot(df, aes_string(x="site.name", y=measure))
  
  # Add those percentile lines
  pg2=pg1 + compile_percentiles(percentiles, measure, qcols)
  #pg2=pg1
  #p <- percentiles[1,]
  #pg2=pg2 + geom_hline(aes_string(yintercept=measure), data=p)
  #p <- percentiles[2,]
  #pg2=pg2 + geom_hline(aes_string(yintercept=measure), data=p)
  #p <- percentiles[3,]
  #pg2=pg2 + geom_hline(aes_string(yintercept=measure), data=p)
  
  
  # Add main plot
  # - violin plot + boxplot for all the data
  # - jitter plot for each site (adjust the color)
  # - x and y labels
  pg3=pg2 + 
    geom_violin(aes(x=global), color="gray50") + 
    geom_boxplot(aes(x=global), width=.1, fill="gray50", outlier.size=0) + 
    geom_jitter(aes(color=site.name), position = position_jitter(width = .1)) + 
    scale_color_manual(values=c(brewer.pal(4,"Dark2"), cmi_site_colors_ramp(nsites))) + 
    ylab(desc) +
    xlab("") 
  
  # Add the y-range limit and the outlier points on the maximum of the range
  # only if there are any outliers
  pg4=pg3
  #if (nrow(lst.outlier.std[[measure]]) > 0) {
  #  pg4=pg4 + 
  #    geom_jitter(aes(x=site.name, y=new), 
  #                data=lst.outlier.std[[measure]], 
  #                position = position_jitter(width = .1), 
  #                shape=8)
  #}
  pg4=pg4 + 
    ylim(df.range.iqr[[measure]])
  
  # Below assumes that you are doing this with a default axis (sites on x, data on y)
  pg5=pg4 + 
    theme_bw() +
    theme(axis.title.x      = element_text(family = "Times", face = "plain", size=12)) +
    theme(axis.title.y      = element_text(family = "Times", face = "plain", size=12, angle=90)) +
    theme(axis.text.x       = element_text(family = "Times", face = "plain", size=10, vjust=0.5, angle=45)) +
    theme(axis.text.y       = element_text(family = "Times", face = "plain", size=10, angle=90)) +
    theme(axis.ticks.length = unit(.15, "lines")) +
    theme(axis.ticks.margin = unit(.15,"lines")) +
    theme(plot.margin       = unit(c(0.25, 0.25,0.25,0.25), "lines"))+
    theme(legend.position   = "none")
  pg=pg5
  print(pg)
  
  #     # Let's now redo the previous plot but with the coordinates flipped
  #     # This also requires changing up the angle of the x-axis
  #     pg5=pg4 + 
  #         theme_bw() +
  #         theme(axis.title.x      = element_text(family = "Times", face = "plain", size=14)) +
  #         theme(axis.title.y      = element_blank()) +
  #         theme(axis.text.x       = element_text(family = "Times", face = "plain", size=12)) +
  #         theme(axis.text.y       = element_text(family = "Times", face = "plain", size=12, angle=0)) +
  #         theme(axis.ticks        = element_blank()) +
  #         theme(plot.margin       = unit(c(0.25, 0.25,0.25,0.25), "lines")) +
  #         theme(legend.position   = "none")
  #     pg6=pg5 + 
  #         coord_flip()
  #     pg=pg6
  #     print(pg)
  
  #ggsave("plot_option02.png", pg, height=2.5, width=5)    
  
  #readline("continue?")
  cat("\n\n\n\n")
}




#' ### KEEP OUTLIERS
#+ plot-qc-yes-outliers, fig.width=15, fig.height=10
for (i in 1:nrow(mdf)) {
  measure <- as.character(mdf$measure[i])
  desc <- as.character(mdf$description[i])
  
  # ### Option 1
  # First, I'll plot ones with 1%, 5%, 25%, and 50% percentile lines
  pg1=ggplot(df, aes_string(x="site.name", y=measure))
  
  # Add those percentile lines
  pg2=pg1 + compile_percentiles(percentiles, measure, qcols)
  #pg2=pg1
  #p <- percentiles[1,]
  #pg2=pg2 + geom_hline(aes_string(yintercept=measure), data=p)
  #p <- percentiles[2,]
  #pg2=pg2 + geom_hline(aes_string(yintercept=measure), data=p)
  #p <- percentiles[3,]
  #pg2=pg2 + geom_hline(aes_string(yintercept=measure), data=p)
  
  
  # Add main plot
  # - violin plot + boxplot for all the data
  # - jitter plot for each site (adjust the color)
  # - x and y labels
  pg3=pg2 + 
    geom_violin(aes(x=global), color="gray50") + 
    geom_boxplot(aes(x=global), width=.1, fill="gray50", outlier.size=0) + 
    geom_jitter(aes(color=site.name), position = position_jitter(width = .1)) + 
    scale_color_manual(values=c(brewer.pal(4,"Dark2"), cmi_site_colors_ramp(nsites))) + 
    ylab(desc) +
    xlab("") 
  
  # Add the y-range limit and the outlier points on the maximum of the range
  # only if there are any outliers
  pg4=pg3
  #if (nrow(lst.outlier.std[[measure]]) > 0) {
  #  pg4=pg4 + 
  #    geom_jitter(aes(x=site.name, y=new), 
  #                data=lst.outlier.std[[measure]], 
  #                position = position_jitter(width = .1), 
  #                shape=8)
  #}
  #pg4=pg4 + 
  #  ylim(df.range.iqr[[measure]])
  
  # Below assumes that you are doing this with a default axis (sites on x, data on y)
  pg5=pg4 + 
    theme_bw() +
    theme(axis.title.x      = element_text(family = "Times", face = "plain", size=12)) +
    theme(axis.title.y      = element_text(family = "Times", face = "plain", size=12, angle=90)) +
    theme(axis.text.x       = element_text(family = "Times", face = "plain", size=10, vjust=0.5, angle=45)) +
    theme(axis.text.y       = element_text(family = "Times", face = "plain", size=10, angle=90)) +
    theme(axis.ticks.length = unit(.15, "lines")) +
    theme(axis.ticks.margin = unit(.15,"lines")) +
    theme(plot.margin       = unit(c(0.25, 0.25,0.25,0.25), "lines"))+
    theme(legend.position   = "none")
  pg=pg5
  print(pg)
  
  #     # Let's now redo the previous plot but with the coordinates flipped
  #     # This also requires changing up the angle of the x-axis
  #     pg5=pg4 + 
  #         theme_bw() +
  #         theme(axis.title.x      = element_text(family = "Times", face = "plain", size=14)) +
  #         theme(axis.title.y      = element_blank()) +
  #         theme(axis.text.x       = element_text(family = "Times", face = "plain", size=12)) +
  #         theme(axis.text.y       = element_text(family = "Times", face = "plain", size=12, angle=0)) +
  #         theme(axis.ticks        = element_blank()) +
  #         theme(plot.margin       = unit(c(0.25, 0.25,0.25,0.25), "lines")) +
  #         theme(legend.position   = "none")
  #     pg6=pg5 + 
  #         coord_flip()
  #     pg=pg6
  #     print(pg)
  
  #ggsave("plot_option02.png", pg, height=2.5, width=5)    
  
  #readline("continue?")
  cat("\n\n\n\n")
}


