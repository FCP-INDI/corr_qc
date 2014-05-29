#!/usr/bin/env Rscript

# This script reads in motion outputs from CPAC and gets it into my framework

library(plyr)

#setwd("~/zarrar")
#setwd("~/Dropbox/Research/cmi/")
#setwd("qc/corr.scripts")

#script.dir <- dirname(sys.frame(1)$ofile)
#setwd(script_dir)

source("10_import_raw_spreadsheet.R")

qc.measures <- c("MeanFD", "PercentFD_greater_than_0.20")


#-------------
# READ IN DATA
#-------------

# Paths
in_qc_file  <- "../corr.qc/raw/epi_temporal/nsd_motion_4669.csv"
#infiles     <- Sys.glob("../corr.qc/raw/epi_temporal/*_params.csv")
in_raw_file <- "../corr.qc/allsites_rest_qc_ips.csv"
in_exc_file <- "../corr.qc/allsites_rest_qc_exclude.txt"

# Read in QC measures
#df1         <- ldply(infiles, read.csv, row.names=NULL)
## Fix the columns because they were saved improperly
#df2         <- df1[,-ncol(df1)]
#colnames(df2) <- colnames(df1)[-1]
qc1         <- read.csv(in_qc_file)
qc2         <- qc1

# Read in complete list of subjects
info        <- read_all_rest(in_raw_file, in_exc_file)


#-------------
# Filter Data
#-------------

# Keep only those columns that you want
qc3         <- subset(qc2, select=c("Subject", "Scan", qc.measures))

# Remove subjects with Exclude in the ID
# and remove subjects starting with sub
bad_ids1    <- grepl("_Exclude_", as.character(qc3$Subject))
bad_ids2    <- grepl("^sub[0-9]*", as.character(qc3$Subject))
qc4         <- qc3[!bad_ids1 & !bad_ids2,]

# Split subject column into subject and session
tmp         <- strsplit(as.character(qc4$Subject), "_")
all(sapply(tmp, length)==3) # double check
subject     <- sapply(tmp, function(x) x[1])
session     <- sapply(tmp, function(x) x[3])
qc5         <- cbind(subid=subject, session=session, qc4[,-1])

# Extract the scan id
scan        <- as.character(qc5$Scan)
scan        <- sub("rest_", "", scan)
scan        <- sub("_rest", "", scan)
scan        <- sub("_1", "", scan)  # something weird is happening here?
scan        <- as.numeric(scan)
qc6         <- cbind(qc5[,1:2], scan=scan, qc5[,4:ncol(qc5)])

# Ensure qc columns are numeric
qc7         <- qc6
for (measure in qc.measures) {
  qc7[[measure]] <- as.numeric(as.character(qc7[[measure]]))
}


#-----------------------------
# INCORPORATE SITE INFORMATION
#-----------------------------

# Merge the QC matrix with site information
df1         <- merge(info, qc7, by=c("subid", "session", "scan"), all.x=TRUE)

# Restrict to only those scans where the raw data exists and is preprocessed
df2         <- df1[df1$raw & df1$preprocessed, ]

# Check if any derivatives are missing
# and exclude (for now BUT WARN)
all_bad_inds <- sapply(qc.measures, function(measure) {
  x         <- as.character(df2[[measure]])
  bad_inds  <- is.na(x) | x == ""
  if (any(bad_inds)) cat(sprintf("%s has %i missing values\n", measure, sum(bad_inds)))
  return(bad_inds)
})
keep_inds   <- rowSums(all_bad_inds) == 0 # exclude any subject with missing QC measure
if (any(!keep_inds)) {
  cat("...removing", sum(!keep_inds), "subjects due to missing data", "\n")
  write.csv(df2[!keep_inds,], file="../corr.qc/bad_qc_filt_motion.csv")
}
df3         <- df2[keep_inds,]

# TODO: convert the site labels to abbreviated version!
# TODO: add site, site.name (where site => site.name and then site => 1:nsites)

# Only include the subset of measures of interest
df4         <- subset(df3, select=c("site", "site.name", "uniqueid", "subid", "session", "scan", qc.measures))

# Ensure qc columns are numeric
df5         <- df4
for (measure in qc.measures) {
  df5[[measure]] <- as.numeric(as.character(df4[[measure]]))
}

# Adjust site.name and number for NKI scans
nki_inds    <- df5$site.name == "NKI"
inds_s1     <- df5$scan[nki_inds] == 645
inds_s2     <- df5$scan[nki_inds] == 1400
inds_s3     <- df5$scan[nki_inds] == 2500
site_names  <- as.character(df5$site.name)
site_names[nki_inds][inds_s1] <- "NKI 1"
site_names[nki_inds][inds_s2] <- "NKI 2"
site_names[nki_inds][inds_s3] <- "NKI 3"
df5$site.name<- factor(site_names)

# Add on a 'global' site column, for any figure combining all the sites together
df5$global   <- factor(rep("All Sites", nrow(df5)))


#-----
# SAVE
#-----

# Save
df           <- df5
write.csv(df, file="../corr.qc/qc_filt_epi_motion.csv")

