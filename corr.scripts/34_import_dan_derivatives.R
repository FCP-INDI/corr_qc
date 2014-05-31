#!/usr/bin/env Rscript

library(plyr)

#setwd("~/zarrar")
#setwd("~/Dropbox/Research/cmi/")
#setwd("qc/corr.scripts")

#script.dir <- dirname(sys.frame(1)$ofile)
#setwd(script_dir)


#-------------
# READ IN DATA
#-------------

in_qc_file  <- "../corr.qc/raw/epi_temporal/FINAL_func_qc.csv"
in_site_file<- "../corr.qc/subject_site_info.csv"

# Read in QC measures
qc1         <- read.csv(in_qc_file)

# QC measures to look for
derivatives <- c("falff", "reho", "vmhc")
measures <- c("50", "75", "90", "fwhm", "mean")
qc.measures <- expand.grid(d=derivatives, m=measures) 
qc.measures <- paste(qc.measures[,1], qc.measures[,2], sep="_")
#qc.measures <- colnames(qc1)[-1]


#-------------
# Filter Data
#-------------

cat("Filtering QC Data\n")

# Split subject column into subject and session
tmp         <- strsplit(as.character(qc1$Subject), "_")
all(sapply(tmp, length)==3) # double check
subject     <- sapply(tmp, function(x) x[1])
session     <- sapply(tmp, function(x) x[3])
qc2         <- cbind(subid=subject, session=session, qc1[,!("subject" == colnames(qc1))])

# Extract the scan id
scan        <- as.character(qc2$Scan)
scan        <- sub("rest_", "", scan)
scan        <- sub("_rest", "", scan)
scan        <- sub("_1", "", scan)  # something weird is happening here?
scan        <- as.numeric(scan)
qc3         <- cbind(qc2[,1:2], scan=scan, qc2[,3:ncol(qc2)])

# ?Should I save the original subject field as unique_id?

# Ensure qc columns are numeric
qc4         <- qc3
for (measure in qc.measures) {
  qc4[[measure]] <- as.numeric(as.character(qc3[[measure]]))
}

# Remove any subjects with missing data
bad_inds    <- apply(subset(qc4, select=qc.measures), 2, function(x) {
  x         <- as.character(x)
  bad_inds  <- is.na(x) | x == ""
  return(bad_inds)
})
keep_inds   <- rowSums(bad_inds) == 0
qc5         <- qc4[keep_inds,]
nbad        <- sum(rowSums(bad_inds) > 0)
cat(sprintf("...removed %i subjects with missing data\n", nbad))
all(!sapply(qc.measures, function(m) any(is.na(qc5[m])))) ## double check


#-----------------------------
# INCORPORATE SITE INFORMATION
#-----------------------------

# Load in the site information
info    <- read_all_rest(in_raw_file, in_exc_file)

# Merge the QC matrix with site information
df1      <- merge(info, qc5, by=c("subid", "session", "scan"), all.x=TRUE)

# Restrict to only those scans where the raw data exists and is preprocessed
df2         <- df1[df1$raw & df1$preprocessed, ]

# Check if any measures are missing
# and exclude (for now BUT WARN)
all_bad_inds <- sapply(qc.measures, function(measure) {
  x         <- as.character(df2[[measure]])
  bad_inds  <- is.na(x) | x == ""
  if (any(bad_inds)) cat(sprintf("%s has %i missing values\n", measure, sum(bad_inds)))
  return(bad_inds)
})
keep_inds   <- rowSums(all_bad_inds) == 0 # exclude any subject with missing QC measure
cat("...removing", sum(!keep_inds), "subjects due to missing data", "\n")
df3         <- df2[keep_inds,]

# Only include the subset of measures of interest
df4         <- subset(df3, select=c("site", "site.name", "uniqueid", "subid", "session", "scan", qc.measures))

# Ensure qc columns are numeric (again...)
df5         <- df4
for (measure in qc.measures) {
  df5[[measure]] <- as.numeric(as.character(df4[[measure]]))
}

# Adjust site.name and number for NKI scans
df5         <- add_nki_samples(df5)

# Add on a 'global' site column, for any figure combining all the sites together
df5$global   <- factor(rep("All Sites", nrow(df5)))


#-----
# SAVE
#-----

df <- df5
write.csv(df, file="../corr.qc/qc_filt_epi_derivatives.csv")

