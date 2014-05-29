#!/usr/bin/env Rscript

# This script gathers the subjects, sessions, and scans for the CoRR data

library(plyr)

pipedir <- "/data/Projects/CoRR/preproc/output/pipeline_corr_qc_preproc"
ids     <- list.files(pattern="[0-9]*_session_*", pipedir)

# check that have usable anatomical scans
# remove subjects that have no anatomical
usable.anat <- laply(ids, function(id) {
    scandirs<- list.files(file.path(pipedir, id, "anatomical_reorient"))
    return( length(scandirs) > 0 )
}, .progress="text")
ids     <- ids[usable.anat] ## all are usable

# check that have usable functional scans
# remove subjects that have no functional
usable.func <- laply(ids, function(id) {
    scandirs<- list.files(file.path(pipedir, id, "mean_functional"))
    return( length(scandirs) > 0 )
}, .progress="text")
ids     <- ids[usable.func]

# add scans
ids2    <- llply(ids, function(id) {
    scandirs<- list.files(file.path(pipedir, id, "mean_functional"))
    scans   <- sub("_rest", "", sub("_scan_rest_", "", scandirs))
    id2     <- paste(id, "scan", scans, sep="_")
    return(id2)
}, .progress="text")
ids2    <- unlist(ids2)

# remove any scan with Exclude in the name or starting with sub
bad_inds1   <- grepl("_Exclude_", ids2)
bad_inds2   <- grepl("^sub", ids2)
keep_inds   <- !(bad_inds1 | bad_inds2)
ids2        <- ids2[keep_inds]

# split information
ids3    <- strsplit(ids2, "_")  # list of lists

# TEMPORARY FIX: 645, 1400, 2500 data
err_inds<- which(sapply(ids3, length) == 6)
ids3    <- lapply(ids3, function(x) {
    if (length(x) == 6) {
        x <- x[1:5]
    }
    return(x)
})

# remove any subjects with not the right number of splits
if (any(sapply(ids3, length)!=5)) stop("unknown subject folder id")

# reshape as data-frame
df      <- ldply(ids3, function(id3) {
    c(subject=id3[1], session=id3[3], scan=id3[5])
})
df$session  <- as.integer(df$session)
df$scan     <- as.integer(df$scan)

# this is another ghetto fix. remove the duplicate scan labels
# find nki subjects and remove those with scan 1
rm_subs <- unique(df$subject[err_inds])
inds    <- which(df$subject %in% rm_subs)
rm_inds <- inds[df$scan[inds]==1]
df      <- df[-rm_inds,]

# save
if (!file.exists("corr.qc")) dir.create("corr.qc")
write.csv(df, file="corr.qc/subject_scan_info.csv")

# save all
if (!file.exists("corr.qc/subject_lists")) dir.create("corr.qc/subject_lists")
d_ply(df, .(session, scan), function(sdf) {
    cat("Session", sdf$session[1], "-", "Scan", sdf$scan[1], "\n")
    outfile <- sprintf("corr.qc/subject_lists/session%02i_scan%02i.txt", sdf$session[1], sdf$scan[1])
    write.table(as.character(sdf$subject), file=outfile, row.names=F, col.names=F, quote=F)
})

# save unique sessions and scans
df <- read.csv("corr.qc/subject_scan_info.csv")
ssdf <- unique(df[,3:4]) 
ssdf <- ssdf[order(ssdf$session, ssdf$scan),]
write.csv(ssdf, file="corr.qc/scan_info.csv")
