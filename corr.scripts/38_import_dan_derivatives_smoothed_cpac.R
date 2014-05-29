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

in_qc_file  <- "../corr.qc/raw/epi_derivative/nsd_derivs_cpac_smoothed_4554.csv"
in_site_file<- "../corr.qc/subject_site_info.csv"

# Read in QC measures
df1         <- read.csv(in_qc_file)

# QC measures to look for
## derivatives <- c("falff", "reho", "vmhc")
## measures <- c("50", "75", "90", "fwhm", "mean")
## qc.measures <- expand.grid(d=derivatives, m=measures) 
## qc.measures <- paste(qc.measures[,1], qc.measures[,2], sep="_")
qc.measures <- colnames(df1)[-1]


#-------------
# Filter Data
#-------------

cat("Filtering QC Data\n")

# Remove subjects with Exclude in the ID
# and remove subjects starting with sub
bad_ids1    <- grepl("_Exclude_", as.character(df1$subject))
bad_ids2    <- grepl("^sub[0-9]*", as.character(df1$subject))
df2         <- df1[!bad_ids1 & !bad_ids2,]
cat(sprintf("...removed %i with bad ids\n", sum(bad_ids1 | bad_ids2)))

# Split subject column into subject / session / scan
tmp         <- strsplit(as.character(df2$subject), "_")
tmp_len     <- sapply(tmp, length)
if (any(!(tmp_len %in% 7:8))) stop("weird subject fields")
## regular ones have length 7
##  "0026017" "session" "4"       "scan"    "rest"    "3"       "rest"
## fix ones with length 8 to have length 7
##  "3201815" "session" "1"       "scan"    "rest"    "1400"    "1" "rest"  
tmp[tmp_len==8] <- lapply(tmp[tmp_len==8], function(x) x[-7])
all(sapply(tmp, length)==7) # double check
## finally can do the splitting
subject     <- sapply(tmp, function(x) x[1])
session     <- sapply(tmp, function(x) x[3])
scan        <- sapply(tmp, function(x) x[6])
df3         <- cbind(subject=subject, session=session, scan=scan, df2[,-1])

# ?Should I save the original subject field as unique_id?

# Ensure qc columns are numeric
df4         <- df3
for (measure in qc.measures) {
  df4[[measure]] <- as.numeric(as.character(df3[[measure]]))
}

# Remove any subjects with missing data
bad_inds    <- apply(subset(df4, select=qc.measures), 2, function(x) {
  x         <- as.character(x)
  bad_inds  <- is.na(x) | x == ""
  return(bad_inds)
})
keep_inds   <- rowSums(bad_inds) == 0
df5         <- df4[keep_inds,]
nbad        <- sum(rowSums(bad_inds) > 0)
cat(sprintf("...removed %i subjects with missing data\n", nbad))
all(!sapply(qc.measures, function(m) any(is.na(df5[m])))) ## double check

## Remove extra NKI subject scans
#nki_subs    <- as.character(df5$subject[df5$scan == 645])
#inds        <- which(df5$subject %in% nki_subs)
#rm_inds     <- inds[df5$scan[inds]==1]
#df7         <- df6[-rm_inds,]
#cat(sprintf("...removed %i nki subjects with duplicate scans\n", length(rm_inds)))

#-----------------------------
# INCORPORATE SITE INFORMATION
#-----------------------------

# Load in the site information
info    <- read.csv(in_site_file, colClasses="character")
info    <- info[,-1]

# Merge the QC matrix with site information
df      <- merge(info, df5, by="subject")

# Adjust site.name and number for NKI scans
nki_inds    <- df$site.name == "NKI"
inds_s1     <- df$scan[nki_inds] == 645
inds_s2     <- df$scan[nki_inds] == 1400
inds_s3     <- df$scan[nki_inds] == 2500
site_names  <- as.character(df$site.name)
site_names[nki_inds][inds_s1] <- "NKI 1"
site_names[nki_inds][inds_s2] <- "NKI 2"
site_names[nki_inds][inds_s3] <- "NKI 3"
df$site.name<- factor(site_names)

# Add on a 'global' site column, for any figure combining all the sites together
df$global   <- factor(rep("All Sites", nrow(df)))


#-----
# SAVE
#-----

write.csv(df, file="../corr.qc/qc_filt_epi_derivatives_smoothed_cpac.csv")

