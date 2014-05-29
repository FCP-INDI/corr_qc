#' ## Outliers
#' Sometimes extreme data-points can skew the plot 
#' and make it difficult to see the spread of the data.
#' 
#' Here I have two functions that use either standard-deviation
#' or IQR to determine if a point is an outlier (6 std or 6 IQR).
#' 
#' I will avoid plotting those outlier points by setting the axis
#' to only the points that I want. ggplot will complain but let her (poor baby).
#' 
#' Outliers will be plotted as astericks at the maximum point in the limited range.
#+ plot-find-outliers
# functions
outlier.std <- function(x) scale(x)
outlier.iqr <- function(x) (x - median(x))/IQR(x)
# outlier vals
df.outlier.std <- ddply(df, .(global), colwise(outlier.std, qc.measures))
df.outlier.iqr <- ddply(df, .(global), colwise(outlier.iqr, qc.measures))
# new ranges
df.range.std <- as.data.frame(sapply(qc.measures, function(m) 
  range(df[[m]][abs(df.outlier.std[[m]]) < 6])) * 1.1
)
df.range.iqr <- as.data.frame(sapply(qc.measures, function(m) 
  range(df[[m]][abs(df.outlier.iqr[[m]]) < 6])) * 1.1
)
# plot outliers
lst.outlier.std <- lapply(qc.measures, function(m) {
  pos.df <- df[df.outlier.std[[m]]>6,]
  neg.df <- df[df.outlier.std[[m]]<(-6),]
  
  odf <- rbind(pos.df, neg.df)
  odf <- subset(odf, select=c("subject", "site", "site.name", "session", "scan"))
  odf$original <- c(pos.df[[m]], neg.df[[m]])
  odf$new <- rep(rev(df.range.std[[m]]), c(nrow(pos.df), nrow(neg.df)))
  odf$new <- odf$new * 0.99
  
  return(odf)
})
lst.outlier.iqr <- lapply(qc.measures, function(m) {
  pos.df <- df[df.outlier.iqr[[m]]>6,]
  neg.df <- df[df.outlier.iqr[[m]]<(-6),]
  
  odf <- rbind(pos.df, neg.df)
  odf <- subset(odf, select=c("subject", "site", "site.name", "session", "scan"))
  odf$original <- c(pos.df[[m]], neg.df[[m]])
  odf$new <- rep(rev(df.range.iqr[[m]]), c(nrow(pos.df), nrow(neg.df)))
  odf$new <- odf$new * 0.99
  
  return(odf)
})
names(lst.outlier.std) <- qc.measures
names(lst.outlier.iqr) <- qc.measures
