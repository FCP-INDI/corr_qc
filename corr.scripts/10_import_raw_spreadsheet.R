# This doesn't exactly import the raw spreadsheet, simply offers functions to do so

read_all_anat <- function(df_file, exclude_file) {
  df <- read.csv(df_file)
  exclude_subs <- as.character(read.table(exclude_file)[,])
  format_all_anat_df(df, exclude_subs)
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

read_all_rest <- function(df_file, exclude_file) {
  df <- read.csv(df_file)
  exclude_subs <- as.character(read.table(exclude_file)[,])
  format_all_rest_df(df, exclude_subs)
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
