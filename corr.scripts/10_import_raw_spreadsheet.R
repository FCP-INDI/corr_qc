# This doesn't exactly import the raw spreadsheet, simply offers functions to do so

add_site_names <- function(df, site_info) {
    info    <- site_info
    isubjs  <- as.numeric(as.character(df$subid))
    info    <- info[as.character(info$ID.Range) != "",]

    #' I can loop through each subject id and see which ID.Range it is in (i.e., the row)
    #' I might want to represent the ID Range differently?
    #' let's then get the min and max
    #+ min-max
    id.ranges <- gsub(",", "", as.character(info$ID.Range))
    id.ranges <- strsplit(id.ranges, "-")
    id.ranges <- t(sapply(id.ranges, as.integer))
    colnames(id.ranges) <- c("start", "end")

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
    site.names  <- as.character(info$CoRR.Plot.Label[isites])
    usite.names <- unique(site.names)

    head(isites)
    head(site.names)

    #df$site.org <- df$site
    df$site <- isites
    df$site.name <- site.names
    
    return(df)
}

add_nki_samples <- function(df) {
    nki_inds    <- df$site.name == "NKI"
    inds_s1     <- df$scan[nki_inds] == 645
    inds_s2     <- df$scan[nki_inds] == 1400
    inds_s3     <- df$scan[nki_inds] == 2500
    site_names  <- as.character(df$site.name)
    site_names[nki_inds][inds_s1] <- "NKI 1"
    site_names[nki_inds][inds_s2] <- "NKI 2"
    site_names[nki_inds][inds_s3] <- "NKI 3"
    df$site.name<- factor(site_names)
    return(df)
}

read_all_anat <- function(df_file, exclude_file, site_info_file="../corr.qc/indi_ids_by_site.csv") {
  df <- read.csv(df_file)
  exclude_subs <- as.character(read.table(exclude_file)[,])
  df <- format_all_anat_df(df, exclude_subs)
  
  site_info <- read.csv(site_info_file)
  df <- add_site_names(df, site_info)
  
  return(df)
}

format_all_anat_df <- function(df, exclude_subs=NULL) {
  # Lower-case column names for easier interaction
  colnames(df) <- tolower(colnames(df))
  
  # Rename columns again for easier use
  old_cols <- c("image", "raw_present", "symlink", "anat_reorient_present", "anat_mni_present")
  new_cols <- c("scan", "raw", "symlink", "reorient", "mni")
  if (!all(old_cols %in% colnames(df))) stop("mismatch in column names")
  inds     <- sapply(old_cols, function(col) which(colnames(df) == col))
  colnames(df)[inds] <- new_cols
  
  # Remove characters from session column
  session <- sub("session_", "", as.character(df$session))
  session <- as.numeric(session)
  if (any(is.na(session))) {
    print(df[is.na(session),])
    stop("error in some session ids")
  }
  df$session <- session
  
  # Remove characters from image column (also rename to scan)
  scan <- sub("anat_", "", as.character(df$scan))
  scan <- as.numeric(scan)
  if (any(is.na(scan))) {
    print(df[is.na(scan),])
    stop("error in some scan/image ids")
  }
  df$scan <- scan
  
  # Relabel Yes/No => True/False
  cols_relabel <- which(sapply(df[1,], function(x) x %in% c("Yes","No")))
  for (col in cols_relabel) {
    vals0 <- as.character(df[[col]])
    vals1 <- vals0 == "Yes"
    df[[col]] <- vals1
  }
  
  # Exclude subjects
  if (!is.null(exclude_subs)) {
    bad_inds <- df$uniqueid %in% exclude_subs
    if (sum(bad_inds) != length(exclude_subs)) {
      cat("found %i subjects to remove but expected to remove %i subjects\n", sum(bad_inds), length(exclude_subs))
    }
    df <- df[!bad_inds,]
  }
  
  return(df)
}

read_all_rest <- function(df_file, exclude_file, site_info_file="../corr.qc/indi_ids_by_site.csv") {
  df <- read.csv(df_file)
  exclude_subs <- as.character(read.table(exclude_file)[,])
  df <- format_all_rest_df(df, exclude_subs)
  
  site_info <- read.csv(site_info_file)
  df <- add_site_names(df, site_info)
  
  return(df)
}

format_all_rest_df <- function(df, exclude_subs=NULL) {
  # Lower-case column names for easier interaction
  colnames(df) <- tolower(colnames(df))
  
  # Rename columns again for easier use
  old_cols <- c("image", "raw_present", "rest_functional_brain_mask_present", "preprocessed_present", "functional_mni_present", "raw_reho.present", "raw_vmhc_present", "raw_falff_present")
  new_cols <- c("scan", "raw", "brain_mask", "preprocessed", "mni", "reho", "vmhc", "falff")
  if (!all(old_cols %in% colnames(df))) stop("mismatch in column names")
  inds     <- sapply(old_cols, function(col) which(colnames(df) == col))
  colnames(df)[inds] <- new_cols
  
  # Remove characters from session column
  session <- sub("session_", "", as.character(df$session))
  session <- as.numeric(session)
  if (any(is.na(session))) {
    print(df[is.na(session),])
    stop("error in some session ids")
  }
  df$session <- session
  
  # Remove characters from image column (also rename to scan)
  scan <- sub("rest_", "", as.character(df$scan))
  scan <- sub("_1", "", scan)
  scan <- as.numeric(scan)
  if (any(is.na(scan))) {
    print(df[is.na(scan),])
    stop("error in some scan/image ids")
  }
  df$scan <- scan
  
  # Relabel Yes/No => True/False
  cols_relabel <- which(sapply(df[1,], function(x) x %in% c("Yes","No")))
  for (col in cols_relabel) {
    vals0 <- as.character(df[[col]])
    vals1 <- vals0 == "Yes"
    df[[col]] <- vals1
  }
  
  # Exclude subjects
  if (!is.null(exclude_subs)) {
    bad_inds <- df$uniqueid %in% exclude_subs
    if (sum(bad_inds) != length(exclude_subs)) {
      cat("found %i subjects to remove but expected to remove %i subjects\n", sum(bad_inds), length(exclude_subs))
    }
    df <- df[!bad_inds,]
  }
  
  return(df)
}

