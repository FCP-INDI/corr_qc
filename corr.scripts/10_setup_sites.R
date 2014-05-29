#' Goal of this script is to spit out a file with the list of subjects
#' and the corresponding sites. 
#' 
#' It should also spit out the unique site names etc but this isn't done
#' yet. [ADD TO TODO].

#' ## CLI Subject List
#'
#' To get the subject list, can run the following command on the command-line
#'  cd /data/Projects/CoRR/preproc/output/pipeline_corr_qc_preproc
#'  ls -d *_session_* | \
#'  grep -v EXCLUDE | \
#'  grep -v '^sub' | \
#'  sed -e s/_session_[0-9]*//g -e s/^sub//g > ~/projects/qc/corr.qc/subjects.txt
#'  cd -
#' CONFIRM THIS

#' ## Read in the data
#' 
#' Read in the spatial epi and temporal measures
#' as well as information about subjects and sites for plotting
#+ read
# Hide this baby here
#setwd("/home2/zarrar/projects/qc")
setwd("~/zarrar/qc")
osubjs  <- as.character(read.table("corr.qc/subjects.txt", colClasses="character")[,])
isubjs  <- as.numeric(osubjs)
info    <- read.csv("corr.qc/indi_ids_by_site.csv")
info    <- info[as.character(info$ID.Range) != "",]

#' I can loop through each subject id and see which ID.Range it is in (i.e., the row)
#' I might want to represent the ID Range differently?
#' let's then get the min and max
#+ min-max
id.ranges <- gsub(",", "", as.character(info$ID.Range))
id.ranges <- strsplit(id.ranges, "-")
id.ranges <- t(sapply(id.ranges, as.integer))
colnames(id.ranges) <- c("start", "end")

#' The NKI-RS IDs are a bit off so we will need to extend the IDs all the way
#' to 9,639,203 even though there aren't that many actual IDs
#+ fix-nki
id.ranges[nrow(id.ranges),2] <- 9639203

#' Now we create a vector that is as long as the maximum of the ID ranges
#' Then we assign each value in the vector to the row in the id.ranges,
#' this row number is basically our site ID.
vec.sites <- vector("integer", length=max(id.ranges))
for (i in 1:nrow(id.ranges)) {
  inds    <- id.ranges[i,1]:id.ranges[i,2]
  vec.sites[inds] <- i
}

#' Now what we want is to give each subject a site
#' note that these site numbers correspond to the row
#' in id.ranges and info.
isites    <- vec.sites[isubjs]

#' Weird that only Berlin Mind and Brain Institute isn't inside CORR
#' Should confirm that these ID-Range and Subject IDs match.
info[unique(isites),]

#' Site Names
#' This is temporary! Not sure now I am using a number.
site.names  <- as.character(info[isites,7])
usite.names <- unique(site.names)

#' Now I want to get the names to display for each site
#' Since this is a little confusing and the names are a bit long, I'll do everything by Site #
#+ site-names
# TODO
df        <- data.frame(subject=osubjs, site=isites, site.name=site.names)
head(df)

#' Some of the subjects will repeat in this list. Keep only unique ones.
#+ unique-subjects
df        <- df[!duplicated(df$subject),]

write.csv(df, file="corr.qc/subject_site_info.csv")
#write.csv(#, file="corr.qc/site_info.csv")
