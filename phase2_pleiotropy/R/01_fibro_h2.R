# Phase 2, step 1: reproduce the fibromyalgia SNP heritability with GenomicSEM.

# The paper reports an observed-scale SNP heritability of 10.4% (95% CI 9.8-11)
# from stratified LD score regression with a 52 annotation baseline model. The univariate LDSC
# used here (single LD score) is a different model

source("phase2_pleiotropy/R/helpers.R")

gwas_in   <- "phase2_pleiotropy/data/fibro/GCST90838603.tsv"   # European primary meta analysis
ref_dir   <- "phase2_pleiotropy/data/ref"
munge_dir <- "phase2_pleiotropy/data/munged"
res_dir   <- "phase2_pleiotropy/results"
for (d in c(ref_dir, munge_dir, res_dir, dirname(gwas_in))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(gwas_in)) stop("Copy GCST90838603.tsv into ", dirname(gwas_in), " first.")

# LD scores and the HapMap3 SNP list
find_hm3 <- function() list.files(ref_dir, pattern = "^w_hm3.*snplist$", recursive = TRUE, full.names = TRUE)
if (length(find_hm3()) == 0) {
  tgz <- file.path(ref_dir, "eur_w_ld_chr.tgz")
  options(timeout = 3600)
  download.file("https://zenodo.org/records/14993076/files/eur_w_ld_chr.tgz", tgz, mode = "wb")
  untar(tgz, exdir = ref_dir)
}
hm3_file <- find_hm3()
if (length(hm3_file) == 0) stop("Could not find w_hm3.snplist after extracting the LD score archive.")
hm3_file <- hm3_file[1]
ld_dir <- paste0(dirname(normalizePath(hm3_file)), "/")
n_ld <- length(list.files(ld_dir, pattern = "l2\\.ldscore\\.gz$"))
message("HapMap3 list: ", hm3_file, "\nLD score folder: ", ld_dir, " (", n_ld, " chromosome files; expect 22)")
hm3_snps <- fread(hm3_file)$SNP

# Prepare the GWAS 
prepared <- file.path(munge_dir, "FM_EUR_prepared.txt")
if (!file.exists(prepared)) {
  prepare_sumstats(
    gwas_in, prepared,
    cols = c(SNP = "rsid", A1 = "effect_allele", A2 = "other_allele", BETA = "beta",
             P = "p_value", N = "n_eff", FRQ = "effect_allele_frequency"),
    hm3_snps = hm3_snps)
}

# munge and LD score regression 
library(GenomicSEM)
prepared_abs <- normalizePath(prepared)
hm3_abs <- normalizePath(hm3_file)
old_wd <- setwd(munge_dir)            # munge writes its output into the working directory
tryCatch(
  munge(files = prepared_abs, hm3 = hm3_abs, trait.names = "FM_EUR"),
  finally = setwd(old_wd))

sumstats <- normalizePath(file.path(munge_dir, "FM_EUR.sumstats.gz"))
ld_out <- ldsc(traits = sumstats, sample.prev = NA, population.prev = NA,
               ld = ld_dir, wld = ld_dir, trait.names = "FM_EUR")

h2  <- ld_out$S[1, 1]
se  <- sqrt(ld_out$V[1, 1])
icp <- tryCatch(ld_out$I[1, 1], error = function(e) NA_real_)
res <- data.frame(trait = "FM_EUR", h2_observed = h2, se = se,
                  ci_low = h2 - 1.96 * se, ci_high = h2 + 1.96 * se, ldsc_intercept = icp)
write.csv(res, file.path(res_dir, "fibro_h2.csv"), row.names = FALSE)
print(res)
cat(sprintf("\nPaper (stratified LDSC, baseline model): 10.4%% (95%% CI 9.8-11.0).\nThis run: %.1f%% (95%% CI %.1f-%.1f).\n",
            100 * h2, 100 * (h2 - 1.96 * se), 100 * (h2 + 1.96 * se)))
