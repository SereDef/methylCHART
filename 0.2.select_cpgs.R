use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.5'
.libPaths(use_library)


# === INFO =====================================================================
# Find overlapping CpGs between 450k and EPICv1 and EPICv2 
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "minfi",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",   # 450K   (hg19)
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",  # EPICv1 (hg19)
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38"    # EPICv2 (hg38)
), force = TRUE)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)


# ---- Load annotations -------------------------------------------------------
anno450k   <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annoEPICv1 <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annoEPICv2 <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# ---- Get CpG IDs ------------------------------------------------------------
cpg450k   <- rownames(anno450k)
cpgEPICv1 <- rownames(annoEPICv1)

# head(rownames(annoEPICv2))

# EPIC v2 probe IDs carry a suffix (e.g. "cg25383568_TC11") because the array
# contains replicate probes targeting the same CpG. 
# Strip the suffix to recover the base cgID (comparable to 450K/EPIv1) & de-duplicate.
cpgEPICv2 <- unique(sub("_.*$", "", rownames(annoEPICv2)))

# # Explore overlap 
# arrays <- list(`450K` = cpg450k, EPICv1 = cpgEPICv1, EPICv2 = cpgEPICv2)
# all_cpgs <- Reduce(union, arrays)
# 
# # Logical membership matrix: rows = CpGs, cols = arrays
# membership <- sapply(arrays, function(x) all_cpgs %in% x)
# 
# # Readable label for which arrays each CpG belongs to, then count
# combo <- apply(membership, 1, function(r) paste(names(arrays)[r], collapse = " & "))
# counts <- sort(table(combo), decreasing = TRUE)
# counts

# 450K & EPICv1 & EPICv2 : 369,543
# EPICv1 & EPICv2 :        351,835
#   450K & EPICv1 :         82,910
#   450K & EPICv2 :         24,484
# EPICv2 only :            177,874    
# EPICv1 only :             61,571
#   450K only :              8,575

# ---- Select CpGs common to all three arrays ----------------------------------
overlapping_cpgs <- Reduce(intersect, list(cpg450k, cpgEPICv1, cpgEPICv2))
length(overlapping_cpgs)

# ---- Keep only true CpG probes -----------------------------------------------
# Drop any SNP (rs), or non-CpG (ch.) probes that survive the intersection

table(substr(overlapping_cpgs, 1, 2))

final_cpgs <- grep("^cg", overlapping_cpgs, value = TRUE)
length(final_cpgs)
# 366,764

writeLines(final_cpgs, "data/cpg_list.txt")

# ---- (optional) drop sex-chromosome probes -----------------------------------
# Note: using the 450K (hg19) annotation (chrX/chrY membership should be the same
#       regardless of genome build and array) 

# autosomal_cpgs <- overlapping_cpgs[
#   !(anno450k[overlapping_cpgs, "chr"] %in% c("chrX", "chrY"))
# ]
# length(autosomal_cpgs)

# ---- (optional) Remove cross-reactive probes ---------------------------------
# Note: maxprobes only ships "450K" and "EPICv1" lists

# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("markgene/maxprobes")
# 
# xreactive <- union(
#   unlist(maxprobes::xreactive_probes(array_type = "450K")),
#   unlist(maxprobes::xreactive_probes(array_type = "EPIC"))
# )

# clean_cpgs <- setdiff(autosomal_cpgs, xreactive)
# length(clean_cpgs)

# ==============================================================================
# ==============================================================================
# # Negative control probes
# mani450k <- minfi::getManifest(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# maniEPIC <- minfi::getManifest(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# # get address of control probes
# all_control <- minfi::getProbeInfo(mani450k, type='Control')
# table(all_control$Type)
# # get only negative controls
# negative_control <- minfi::getControlAddress(mani450k, controlType='NEGATIVE')
