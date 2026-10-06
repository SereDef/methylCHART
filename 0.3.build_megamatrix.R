use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.5'
.libPaths(use_library)

library(bigmemory)
# library(foreach)

options(bigmemory.allow.dimnames = TRUE)

# === INFO =====================================================================
# Save both cohort specific and mega matrix

data_dir <- '~/methylCHART/data'

pheno <- read.csv(file.path(data_dir, 'pheno_data.csv'))
cpg_list <- readLines(file.path(data_dir, 'cpg_list.txt'))

n_probes <- length(cpg_list)
n_samples <- nrow(pheno)

genr_samples <- pheno[pheno$Cohort == 'GENR',   'Sample_ID'] 
alsp_samples <- pheno[pheno$Cohort == 'ALSPAC', 'Sample_ID'] 

stopifnot(length(genr_samples) + length(alsp_samples) == n_samples)

# --- helper: empty file-backed matrix, rows = samples, cols = CpGs ------------
init_matrix <- function(desc, col_names, row_names, loc = data_dir) {
  desc_file <- file.path(loc, paste0(desc,".desc"))
  
  if (file.exists(desc_file)) {
    warning(' * "', desc_file, '" already exists - replacing it.')
    file.remove(desc_file)
    file.remove(gsub('.desc', '.bin', desc_file))
  }
  
  cli::cli_alert_info('Generate methylation matrix: {length(row_names)} x {length(col_names)}')
  
  bigmat <- bigmemory::filebacked.big.matrix(
    nrow = length(row_names),
    ncol = length(col_names),
    type = "double",
    init = NA_real_,
    backingpath = loc,
    descriptorfile = basename(desc_file),
    backingfile = gsub('.desc', '.bin', basename(desc_file)),
    dimnames = list(row_names, col_names)
  )
  
  return(bigmat)
}

# === Generation R =============================================================

genr_dnam_dir <- '~/GENR3/Methylation/Release4'
normalization <- "Functional"

# Initialize file backed matrix
genr_mat <- init_matrix(desc = 'genr_dnam', col_names = cpg_list, row_names = genr_samples)

# Define subsets of data (ie. Period x Array)
genr_dat <- pheno[pheno$Cohort=='GENR',]
subsets <- unique(genr_dat[, c( "Array","Period")])

# helper: find the 
genr_fetch_dnam <- function(period, array) {
  
  cli::cli_rule("Cohort: {.strong Generation R} | Period: {.strong {period}} | Array: {.strong {array}}")
  
  period_dir <- ifelse(period == 'Birth', 'Birth', paste0(period,"y"))
  array_dir  <- ifelse(array == "450k", "450K", array)
  
  if (period == 'Birth' & array == 'EPICv2') {
    data_loc <- file.path(genr_dnam_dir, 'Bridges', normalization)
  } else {
    data_loc <- file.path(genr_dnam_dir, period_dir, array_dir, normalization)
  }
  
  rds_path <- list.files(path = data_loc, pattern = "ALL", full.names = TRUE)
  
  if (length(rds_path) != 1L) {
    cli::cli_alert_warning("No or multiple matching RDS file found, skipping subset.")
    print(rds_path)
    return(NULL)
  }
  
  cli::cli_progress_step("Reading data", spinner = TRUE)
  cpg_data <- readRDS(rds_path)
  cli::cli_progress_done()
  
  return(cpg_data)
}

# foreach::foreach(i = 1:nrow(subsets),
#                  .packages = c("bigmemory", "cli"),
#                  .export = c('subsets', 'genr_mat',
#                              'genr_fetch_dnam', 'genr_sel',
#                              "cpg_list"),
#                  .combine = 'c') %do% { # %dorng% {


for (i in seq_len(nrow(subsets))){

  period <- subsets$Period[i]
  array <- subsets$Array[i]
  
  cpg_data <- genr_fetch_dnam(period, array)
  if (is.null(cpg_data)) next 
  
  want_samples <- genr_dat[genr_dat$Array == array & 
                           genr_dat$Period == period, 'Sample_ID']
  
  # ---- restrict to wanted CpGs & samples, in a defined order ----------------
  avail_cpgs    <- intersect(cpg_list,     rownames(cpg_data))
  avail_samples <- intersect(want_samples, colnames(cpg_data))
  
  if (!length(avail_cpgs) || !length(avail_samples)) {
    cli::cli_alert_warning('Nothing to write (CpGs: {length(avail_cpgs)}, samples: {length(avail_samples)}).')
    rm(cpg_data); next
  }
  
  if (length(avail_cpgs) < n_probes) {
    cli::cli_alert_info('{n_probes - length(avail_cpgs)} of {n_probes} target CpGs missing in this file.')
  }
  if (length(avail_samples) < length(want_samples)) {
    cli::cli_alert_info('{length(want_samples) - length(avail_samples)} expected samples missing in this file.')
  }
    
  # one subset + one transpose -> samples x CpGs, then free space
  cpg_data <- t(cpg_data[avail_cpgs, avail_samples])
  
  cli::cli_alert_success('Writing {nrow(cpg_data)} samples x {ncol(cpg_data)} CpGs.')
  genr_mat[avail_samples, avail_cpgs] <- cpg_data
}
rm(cpg_data, avail_cpgs, avail_samples, genr_samples)

# === ALSPAC ===================================================================
# library(bigmelon) # damn so many dependencies... nah
library(gdsfmt)

aslp_dnam_file <- '~/Alspac_Tempo/methylation/data/betas/common.gds'
# normalization <- "Functional"

# Read pds file structure
aslp_dnam_link <- openfn.gds(aslp_dnam_file)
print(aslp_dnam_link) # See internal tree structure
# |--+ row.names   { Str8 450838, 4.7M }
# |--+ col.names   { Str8 8670, 158.3K }
# \--+ matrix   { Float64 450838x8670, 29.1G }

node <- index.gdsn(aslp_dnam_link, "matrix")

rn = read.gdsn(index.gdsn(aslp_dnam_link, "row.names")) # cpgs
cn = read.gdsn(index.gdsn(aslp_dnam_link, "col.names")) # samples

avail_cpgs    <- intersect(cpg_list,     rn)
avail_samples <- intersect(alsp_samples, cn)

if (!length(avail_cpgs) || !length(avail_samples)) {
  cli::cli_alert_error('Nothing to write (CpGs: {length(avail_cpgs)}, samples: {length(avail_samples)}).')
}

if (length(avail_cpgs) < n_probes) {
  cli::cli_alert_info('{n_probes - length(avail_cpgs)} of {n_probes} target CpGs missing in this file.')
}
if (length(avail_samples) < length(alsp_samples)) {
  cli::cli_alert_info('{length(want_samples) - length(avail_samples)} expected samples missing in this file.')
}

# Initialize file backed matrix
alsp_mat <- init_matrix(desc = 'alsp_dnam', col_names = cpg_list, row_names = alsp_samples)

# Reading all data at once is not memory safe 
# cpg_data <- read.gdsn(index.gdsn(aslp_dnam_link, "matrix")) # Extract data
# rownames(cpg_data) <- rn
# colnames(cpg_data) <- cn
# 
# alsp_mat[avail_samples, avail_cpgs] <- t(cpg_data[avail_cpgs, avail_samples])

# Precompute selections
# Map data file order to target positions
m <- match(rn, cpg_list)   # for each CpG in the file: its position in cpg_list, or NA
row_sel  <- !is.na(m)  
tgt_cols <- m[row_sel]   # target column positions, in file order

s <- match(cn, alsp_samples)
samp_pos <- which(!is.na(s))    # ID positions to read

# ---- Chunk size from a memory budget ----------------------------------------
chunk_size <- 600
chunks     <- split(samp_pos, ceiling(seq_along(samp_pos) / chunk_size))

cli::cli_alert_info('Reading {sum(row_sel)} CpGs x {length(samp_pos)} samples in {length(chunks)} chunk{?s}.')

cli::cli_progress_bar('Writing', total = length(chunks))
for (cols in chunks) {
  col_sel <- logical(length(cn))
  col_sel[cols] <- TRUE
  
  block <- readex.gdsn(node, sel = list(row_sel, col_sel), simplify = "none")  # cpgs x samples
  
  tgt_rows <- match(cn[cols], alsp_samples)
  alsp_mat[tgt_rows, tgt_cols] <- t(block)
  
  rm(block)
  cli::cli_progress_update()
}
cli::cli_progress_done()

if (inherits(alsp_mat, "big.matrix")) bigmemory::flush(alsp_mat)
closefn.gds(aslp_dnam_link)


