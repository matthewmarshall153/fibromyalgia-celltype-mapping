# Shared helpers for Phase 2 
library(data.table)

# Read a GWAS file, rename columns to the names GenomicSEM::munge expects, keep only
# HapMap3 SNPs, and write a clean tab-delimited file.
#   cols: named character vector, our name = column name in the file
#         c(SNP = "rsid", A1 = "effect_allele", A2 = "other_allele", BETA = "beta",
#           P = "p_value", N = "n_eff", FRQ = "effect_allele_frequency")
#   n_value: optional single number used for N when the file has no per variant N
prepare_sumstats <- function(path_in, path_out, cols, hm3_snps, n_value = NULL) {
  need <- unname(cols)
  dt <- fread(path_in, select = need)
  setcolorder(dt, need)
  setnames(dt, need, names(cols))
  n_raw <- nrow(dt)

  dt[, `:=`(BETA = as.numeric(BETA), P = as.numeric(P))]
  dt <- dt[!is.na(SNP) & grepl("^rs[0-9]+$", SNP)]
  dt <- dt[is.finite(BETA) & !is.na(P) & P > 0 & P <= 1]
  dt <- dt[SNP %in% hm3_snps]
  dt[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  dt <- dt[nchar(A1) == 1 & nchar(A2) == 1]          # single base alleles only

  if (!is.null(n_value)) dt[, N := n_value]
  if ("N" %in% names(dt)) setorder(dt, -N) else setorder(dt, P)
  dt <- dt[!duplicated(SNP)]                          # keep the best-powered row per rsID

  if ("FRQ" %in% names(dt)) {
    dt[, MAF := pmin(as.numeric(FRQ), 1 - as.numeric(FRQ))]
    dt[, FRQ := NULL]
  }
  fwrite(dt, path_out, sep = "\t")
  message(sprintf("prepare_sumstats: %s rows read -> %s kept (HapMap3, biallelic SNPs, deduplicated)",
                  format(n_raw, big.mark = ","), format(nrow(dt), big.mark = ",")))
  invisible(list(rows_in = n_raw, rows_out = nrow(dt), path = path_out))
}
