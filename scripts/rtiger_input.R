library(Matrix)
library(vcfR)

#input_dir <- "/home/jool/athaliana_single_cell/workspace/cb_snp_counts/col0xdb1"
#out_dir <-  "/home/jool/athaliana_single_cell/workspace/RTIGER/col0xdb1/samples"
out_dir <- snakemake@output[[1]]
if (!dir.exists(out_dir)) 
  dir.create(out_dir)

sample_IDs <- readLines(snakemake@input[["samples"]])

vcf <- read.vcfR(snakemake@input[["vcf"]])
vcf_fix <- getFIX(vcf)


AD_m <- readMM(snakemake@input[["ad"]])
DP_m <- readMM(snakemake@input[["dp"]])

# parse snp data from base.vcf and convert to zero-based position index
snp_info <- data.frame(
  SeqID = vcf_fix[, "CHROM"],
  Pos   = as.integer(vcf_fix[, "POS"]) - 1,
  RefA  = vcf_fix[, "REF"],
  AltA  = vcf_fix[, "ALT"]
)

REF <- DP_m - AD_m # DP is aggregate REF+ALT; https://cellsnp-lite.readthedocs.io/en/latest/main/manual.html#brief-introduction

stopifnot(length(sample_IDs) == ncol(AD_m))

# Write one tsv per barcode following https://github.com/rfael0cm/RTIGER#preparing-input-data
for (i in seq_along(sample_IDs)) {
  out <- data.frame(
    SeqID = snp_info$SeqID,
    Pos   = snp_info$Pos,
    RefA  = snp_info$RefA,
    RefC  = as.integer(REF[, i]),
    AltA  = snp_info$AltA,
    AltF  = as.integer(AD_m[, i])
  )
  
  write.table(out, file=paste0(sample_IDs[i], ".tsv"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
}