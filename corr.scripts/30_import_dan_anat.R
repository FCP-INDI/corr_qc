#!/usr/bin/env Rscript

#' Currently, imports the spatial anatomical QC measures from Dan into my framework.

#setwd("~/zarrar")
#setwd("~/Dropbox/Research/cmi/")
#setwd("qc/corr.scripts")

#script.dir <- dirname(sys.frame(1)$ofile)
#setwd(script_dir)

source("10_import_raw_spreadsheet.R")

in_qc_file  <- "../corr.qc/raw/anat_spatial/nsd_anat_2618.csv"
in_raw_file <- "../corr.qc/allsites_anat_qc_ips.csv"
in_exc_file <- "../corr.qc/allsites_anat_qc_exclude.txt"

# QC measures to look for
qc.measures <- c("snr", "cnr", "fber", "efc", "qi1", "fwhm")


###
# Read/Format QC Measure DF
###

# Read in QC measures
qc1         <- read.csv(in_qc_file)

# Get the subject/session info
# and exclude any subjects that should be excluded (yeah that's an explanation)
subsplit    <- strsplit(as.character(qc1$subject), "_")
exclude     <- sapply(subsplit, length)!=3
subsplit    <- subsplit[!exclude]
qc1         <- qc1[!exclude,]
subject     <- sapply(subsplit, function(x) x[1])
session     <- sapply(subsplit, function(x) x[3])

# Combine the subject/session info
qc2         <- cbind(subid=subject, session=session, scan=rep(1,nrow(qc1)), qc1[,-1])
qc2         <- qc2[order(qc2$subid, qc2$session),]

# Take only the subset of measures of interest
qc3         <- subset(qc2, select=c("subid", "session", "scan", qc.measures))


###
# Read/Merge/Format All Subject Information
###

# Load in complete subject information
info        <- read_all_anat(in_raw_file, in_exc_file)

# Merge QC DF with Complete Subject DF
# note if something is missing in df3 then NAs will be added in
df1         <- merge(info, qc3, by=c("subid", "session", "scan"), all.x=TRUE)

# Restrict to only those scans where the raw data exists and is not a symlink
df2         <- df1[df1$raw & !df1$symlink, ]

# Check if any derivatives are missing
# and exclude (for now BUT WARN)
 all_bad_inds <- sapply(qc.measures, function(measure) {
  x         <- as.character(df2[[measure]])
  bad_inds  <- is.na(x) | x == ""
  if (any(bad_inds)) cat(sprintf("%s has %i missing values\n", measure, sum(bad_inds)))
  return(bad_inds)
})
keep_inds   <- rowSums(all_bad_inds) == 0 # exclude any subject with missing QC measure
cat("...removing", sum(!keep_inds), "subjects due to missing data", "\n")
write.csv(df2[!keep_inds,], file="../corr.qc/bad_qc_filt_anat.csv")
df3         <- df2[keep_inds,]

# Only include the subset of measures of interest
df4         <- subset(df3, select=c("site", "site.name", "uniqueid", "subid", "session", "scan", qc.measures))

# Ensure qc columns are numeric
df5 <- df4
for (measure in qc.measures) {
  df5[[measure]] <- as.numeric(as.character(df4[[measure]]))
}

# Adjust site.name and number for NKI scans
nki_inds    <- df5$site.name == "NKI"
site_names  <- as.character(df5$site.name)
site_names[nki_inds] <- "NKI 1-3"
df5$site.name<- factor(site_names)

# Add on a 'global' site column, for any figure combining all the sites together
df5$global   <- factor(rep("All Sites", nrow(df5)))

# Save
df <- df5
write.csv(df, file="../corr.qc/qc_filt_anat.csv")
nrow(df)
