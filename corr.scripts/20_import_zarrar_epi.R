#' This script merges the QC measures with the site information
#'
#' TODO: figure out phase encoding for the ghosting

library(plyr)

# Session and Scan Information
#setwd("/home2/zarrar/projects/qc")
#setwd("~/zarrar/qc")
#setwd("~/Dropbox/Research/cmi/qc")

source("10_import_raw_spreadsheet.R")

qc.measures <- c("efc", "fber", "fwhm", "gsr", "dvars", "quality")

in_raw_file <- "../corr.qc/allsites_rest_qc_ips.csv"
in_exc_file <- "../corr.qc/allsites_rest_qc_exclude.txt"

ss_info <- read.csv("../corr.qc/scan_info.csv")[,-1]

# Loop through sessions and scans
# Reading in all the QC information
df <- ddply(ss_info, .(session, scan), function(sdf) {    
  #session <- 1
  #scan    <- 645
  session <- sdf$session
  scan    <- sdf$scan
  
  cat("Session", session, "Scan", scan, "\n")
  
  # Read in the data
  qc_spat <- read.csv(sprintf("../corr.qc/raw/epi_spatial/qc_spatial_epi_session%02i_scan%02i.csv", session, scan), colClasses="character")
  qc_temp <- read.csv(sprintf("../corr.qc/raw/epi_temporal/qc_temporal_epi_session%02i_scan%02i.csv", session, scan), colClasses="character")
  info    <- read_all_rest(in_raw_file, in_exc_file)
  
  # Remove first columns
  qc_spat <- qc_spat[,-1]
  qc_temp <- qc_temp[,-1]
    
  # Fix to subject column
  qc_spat$subject <- sprintf("%07i", as.integer(as.character(qc_spat$subject)))
  qc_temp$subject <- sprintf("%07i", as.integer(as.character(qc_temp$subject)))
  
  # Weird but due to an error on my part
  # Remove any duplicate subject ids
  qc_spat   <- qc_spat[!duplicated(qc_spat$subject),]
  qc_temp   <- qc_temp[!duplicated(qc_temp$subject),]
  
  # Merge the QC matrices
  # Then merge everything else together
  qc      <- merge(qc_spat, qc_temp, by="subject")
  
  # Add the scan to these columns
  qc$session  <- session
  qc$scan     <- scan
  
  df      <- merge(info, qc, by.x=c("subid", "session", "scan"), by.y=c("subject", "session", "scan"))
  #print(head(df, n=2))
  cat(nrow(qc), "=>", nrow(df), "\n")
  if (nrow(qc_spat) != nrow(qc_temp)) {
    cat("error mismatch in number of rows\n")
    cat(nrow(qc_spat), nrow(qc_temp), "\n")
  }
  
  # Let's filter for the columns of interest
  #print(colnames(df))
  use.cols    <- c("uniqueid", "subid", "site", "site.name", "session", "scan", "efc", "fber", "fwhm", "ghost_x", "ghost_y", "dvars", "quality")
  df          <- subset(df, select=use.cols)
  #head(df)
  
  # This will convert the columns that should be numeric to numeric
  exclude.cols <- c("subid", "uniqueid", "site", "site.name")
  include.cols <- colnames(df)[!(colnames(df) %in% exclude.cols)]
  for (col in include.cols) {
    df[[col]] <- as.numeric(as.character(df[[col]]))
  }
    
  # Add on a 'global' site column, for any figure combining all the sites together
  df$global   <- factor(rep("All Sites", nrow(df)))
  
  # Basic site information
  #nsites      <- length(unique(df$site.name))
  #print(unique(df$site))
  
  cat("\n\n")
  
  return(df)
})

# For now, use which side there is more ghosting on average
# to determine phase encoding and which ghost direction to plot
# TODO: double check any that are IS?
df2         <- df
mg          <- ddply(df2, .(site), colwise(mean, .(ghost_x, ghost_y)))
df3         <- ddply(df2, .(site), function(sdf) {
  comp    <- mean(sdf$ghost_x) > mean(sdf$ghost_y)
  if (comp) {
    sdf$gsr <- sdf$ghost_x
  } else {
    sdf$gsr <- sdf$ghost_y
  }
  return(sdf)
}, .progress="text")
df3         <- subset(df3, select=c("uniqueid", "subid", "site", "site.name", "session", "scan", qc.measures, "global"))

# Adjust site.name and number for NKI scans
df4         <- add_nki_samples(df3)

# Done
df          <- df4
write.csv(df, file="../corr.qc/qc_filt_epi.csv")

