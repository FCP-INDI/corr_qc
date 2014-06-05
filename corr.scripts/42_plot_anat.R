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

#setwd("~/zarrar")
#setwd("~/Dropbox/Research/cmi/")
#setwd("qc/corr.scripts")

#script.dir <- dirname(sys.frame(1)$ofile)
#setwd(script_dir)

df <- read.csv("../corr.qc/qc_filt_anat.csv")[,-1]
nsites <- length(unique(df$site))

# Remove IPCAS 4 and IBATRT
df2 <- subset(df, !(site.name %in% c("IPCAS 4", "IBATRT")))
df2$site.name <- factor(df2$site.name)
df <- df2

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

#' Now let's get the percentiles (for all the data)
#+ percentile
percentiles <- apply(subset(df, select=qc.measures), 2, quantile, qvals, na.rm=TRUE)
percentiles <- as.data.frame(cbind(percentiles, qcat, qline, qsize))
percentiles$qline <- as.factor(qline)
percentiles$qcat  <- as.factor(qcat)
print(percentiles)

#' ## Measure Descriptions
#' Associate a detailed description with each measure
#+ setup-measure
mdf <- data.frame(
  measure = qc.measures, 
  description = c(
    "Signal to Noise Ratio (SNR)", 
    "Contrast to Noise Ratio (CNR)", 
    "Foreground to Background Energy Ratio (FBER)", 
    "Entropy Focus Criteria (EFC)", 
    "Percent of Artifact Voxels (Qi1)", 
    "FWHM in mm (FWHM)"
  )
)
print(mdf)
cat("DOUBLE CHECK THAT COLUMNS OF EACH ROW MATCH ABOVE\n")

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

#' Recompute the percentiles for plots without outliers
#+ percentiles-outliers
percentiles.no.outlier <- sapply(qc.measures, function(measure) {
  df[[measure]]
  bad_uids <- as.character(lst.outlier.iqr[[measure]]$uniqueid)
  sub_df <- subset(df, !(uniqueid %in% bad_uids), select=measure, drop=T)
  quantile(sub_df, qvals, na.rm=TRUE)
})
percentiles.no.outlier <- as.data.frame(cbind(percentiles.no.outlier, qcat, qline, qsize))
percentiles.no.outlier$qline <- as.factor(qline)
percentiles.no.outlier$qcat  <- as.factor(qcat)
print(percentiles.no.outlier)


#' ## Visualization of Text
#' A function with all the theme jazz
#+ viz-text
set_themes <- function(family="Times", text.size.x=14, text.size.y=16, title.size=18) {
  family <- "sans"
  pg <- list(
    theme_bw(), 
    theme(axis.title.x      = element_text(family = family, face = "plain", size=title.size)), 
    theme(axis.title.y      = element_text(family = family, face = "plain", size=title.size, angle=90, vjust=0.25)), 
    theme(axis.text.x       = element_text(family = family, face = "plain", size=text.size.x, vjust=0.5, angle=45)), 
    theme(axis.text.y       = element_text(family = family, face = "plain", size=text.size.y, angle=90)), 
    theme(axis.ticks.length = unit(.15, "lines")), 
    theme(axis.ticks.margin = unit(.15,"lines")), 
    theme(plot.margin       = unit(c(0.25, 1, 0.25, 1), "lines")), 
    theme(legend.position   = "none")
  )
  return(pg)
}


#' ## QC Spatial/Anatomical Measures
#' I will be plotting a bunch of different data-sets here. First I will have 
#' all the data, then I will have it remote 3 x the IQR. So two sets of the
#' same plots.
#'
#' ### REMOVE OUTLIERS
#+ plot-qc-no-outliers, fig.width=12, fig.height=8, dpi=100
for (i in 1:nrow(mdf)) {
  measure <- as.character(mdf$measure[i])
  desc <- as.character(mdf$description[i])
  
  # ### Option 1
  # First, I'll plot ones with 1%, 5%, 25%, and 50% percentile lines
  pg1=ggplot(df, aes_string(x="site.name", y=measure))
  
  # Add those percentile lines
  pg2=pg1 + compile_percentiles(percentiles.no.outlier, measure, qcols)  
  
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
  
  # Add the y-range limit
  pg4=pg3
  pg4=pg4 + 
    ylim(df.range.iqr[[measure]])
  
  # Below assumes that you are doing this with a default axis (sites on x, data on y)
  pg5=pg4 + set_themes()
  
  # Plot
  pg=pg5
  print(pg)
  
  #ggsave("plot_option02.png", pg, height=2.5, width=5)    
  
  #readline("continue?")
  cat("\n\n\n\n")
}


#' ### KEEP OUTLIERS
#+ plot-qc-yes-outliers, fig.width=12, fig.height=8, dpi=100
for (i in 1:nrow(mdf)) {
  measure <- as.character(mdf$measure[i])
  desc <- as.character(mdf$description[i])
  
  # ### Option 1
  # First, I'll plot ones with 1%, 5%, 25%, and 50% percentile lines
  pg1=ggplot(df, aes_string(x="site.name", y=measure))
  
  # Add those percentile lines
  pg2=pg1 + compile_percentiles(percentiles, measure, qcols)  
  
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
  #pg4=pg4 + 
  #  ylim(df.range.iqr[[measure]])
  
  # Below assumes that you are doing this with a default axis (sites on x, data on y)
  pg5=pg4 + set_themes()
  
  # Plot
  pg=pg5
  print(pg)
  
  #ggsave("plot_option02.png", pg, height=2.5, width=5)    
  
  #readline("continue?")
  cat("\n\n\n\n")
}




#' ## Correlogram
#' We compare the different anatomical QC measures to try and understand how the different measures might be related. It appears that we have three different types of measures.
#+ plot-cor
library(corrplot)
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))
cmat <- cor(df[,qc.measures])
corrplot(cmat, order="AOE", type="upper", tl.pos="tp", col=col2(200), outline=T)
corrplot(cmat, add=TRUE, type="lower", method="number", order="AOE", col="black", diag=FALSE, tl.pos="n", cl.pos="n")
