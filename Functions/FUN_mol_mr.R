## (1) function to format the gcta format data to TwoSampleMR format -----
gcta_to_2smr = function(gcta_tb, trait_name, type = c("exposure", "outcome"), ...) {
  require(TwoSampleMR)
  gcta_tb |>
    mutate(outcome = trait_name) |>
    format_data(
      type = type,
      effect_allele_col = "A1",
      other_allele_col = "A2",
      phenotype_col = "outcome",
      eaf_col = "freq",
      samplesize_col = "N",
      beta_col = "b",
      se_col = "se",
      pval_col = "p",
      ...
    )
} 

munge_to_2smr = function(munge_tb, trait_name, type = c("exposure", "outcome"), ...) {
  require(TwoSampleMR)
  munge_tb |>
    mutate(outcome = trait_name) |>
    format_data(
      type = type,
      effect_allele_col = "A2",
      other_allele_col = "A1",
      phenotype_col = "outcome",
      eaf_col = "FRQ",
      samplesize_col = "N",
      beta_col = "BETA",
      se_col = "SE",
      pval_col = "P",
      chr_col = "CHR",
      pos_col = "BP",
      ...
    )
} 

## (2) Perform MR on molecular trait ----
## require the exp_dat/out_dat in TwoSampleMR format
## r2 is the r2 threshold for the clumping
## bfile is the ld reference file
## plink_path is the path to plink
mr_mol = function(exp_dat, out_dat, r2 = 0.3, bfile, plink_path, detailed = FALSE) {
  require(tidyverse)
  require(glue)
  require(TwoSampleMR)
  require(MendelianRandomization)
  require(ieugwasr)
  require(gwasglue)

## Harmonize data
  harmdat = harmonise_data(exp_dat, out_dat, action = 2) |>
    mutate(rsid = SNP)

  harmdat = harmdat |> #- remove ambiguous palindromic SNPs
    filter(mr_keep)

# head(harmdat)
# head(exp_dat)
# filter(harmdat, SNP %in% mesnp)
# filter(exp_dat, SNP == mesnp)
# filter(out_dat, SNP %in% mesnp)

## Steiger filtering (no need for lifespan/end-life traits)
  # violate_snp = harmdat |>
  #   steiger_filtering() |>
  #   select(SNP, steiger_dir, steiger_pval) |>
  #   filter(!steiger_dir) |>
  #   pull(SNP)

#- LD clumping
  clumpdat = ld_clump(harmdat, bfile = bfile, plink_bin = plink_path, clump_r2 = r2)
# get snps overlap with reference data.
# check_ref = ieugwasr::ld_matrix(harmdat$SNP, pop = "EUR", bfile = bfile, plink_bin = plink_path) |>
#   colnames() |>
#   str_extract("^rs[0-9]+")
# exp1 = harmdat$SNP[harmdat$SNP %in% check_ref]
# exp2 = harmdat$SNP[!harmdat$SNP %in% check_ref]
# eQTL$SNP |>
#   str_detect("^rs[0-9]+$") |>
#   mean() * All the SNPs here are in "^rs[0-9]+$" format
# dat = harmdat |> # combine the snp not on reference
#   filter(SNP %in% exp2) |>
#   bind_rows(clumpdat)
  dat = clumpdat

# clumpdat = ld_clump(harmdat, bfile = bfile, plink_bin = plink_path, clump_r2 = 0.01)

## Generalized IVW and MR-Egger for the "MendelianRandomization" package --------
  exp = unique(dat$exposure)
  out = unique(dat$outcome)
  dir = directionality_test(dat) |>
    select(snp_r2.exposure:steiger_pval)
  any_violate = any(dat$SNP %in% violate_snp)
  if (nrow(dat) == 1) { ## One instrument: Wald_ratio only
    result = mr_wald_ratio(b_exp = dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure, se_out = dat$se.outcome)
    result = tibble(
      exposure = exp,
      method = "Wald_ratio",
      nsnp = result$nsnp,
      b = result$b,
      se = result$se,
      pval = result$pval,
      intercept = NA,
      intercept_se = NA,
      intercept_pval = NA,
      hetero_Q = NA,
      hetero_pvale = NA,
      egger_b = NA,
      egger_se = NA,
      egger_p = NA
    )
  } else if (nrow(dat) >= 2) { ## more than Two instruments & all have ld info: gIVW
    rho = ieugwasr::ld_matrix(dat$SNP, pop = "EUR", bfile = bfile, plink_bin = plink_path)
    result = MendelianRandomization::mr_ivw(mr_input(
      bx = dat$beta.exposure,
      bxse = dat$se.exposure,
      by = dat$beta.outcome,
      byse = dat$se.outcome,
      cor = rho
    ))
    result = tibble(
      exposure = exp,
      method = "IVW",
      nsnp = result@SNPs,
      b = result@Estimate,
      se = result@StdError,
      pval = result@Pvalue,
      intercept = NA,
      intercept_se = NA,
      intercept_pval = NA,
      hetero_Q = result@Heter.Stat[1],
      hetero_pvale = result@Heter.Stat[2],
      egger_b = NA,
      egger_se = NA,
      egger_p = NA
    )
    if (nrow(dat) > 2) { ## more than Two instruments: MR-Egger
      result2 = MendelianRandomization::mr_egger(mr_input(
        bx = dat$beta.exposure,
        bxse = dat$se.exposure,
        by = dat$beta.outcome,
        byse = dat$se.outcome,
        cor = rho
      ))
      result$intercept = result2@Intercept
      result$intercept_se = result2@StdError.Int
      result$intercept_pval = result2@Pvalue.Int
      result$egger_b = result2@Estimate
      result$egger_se = result2@StdError.Est
      result$egger_p = result2@Pvalue.Est
    }
  }
  result$outcome = out
  result$any_violate_dir_inst = any_violate
  result = cbind(result, dir)
  if (detailed) {
    return(list(result, dat))
  } else {
    return(result)
  }
}

library(purrr)
mr_mol_poss = possibly(mr_mol, otherwise = NULL, quiet = TRUE)

## Function to run GCTA
# run_gcta

# gcta_path = "/home/kying/gcta_1.93.2beta/gcta64"
# bfile = "LD/1000G.EUR.QC_bfile"
# plink_path = "/programs/plink"

# require(glue)
# require(tidyverse)


# harmdat = harmonise_data(fmt_exp, fmt_out, action = 2) |>
#   mutate(rsid = SNP)

# harmdat = harmdat |> #- remove ambiguous palindromic SNPs
#   filter(mr_keep)

# dt = harmdat |>
#   select(
#     SNP,
#     A1 = effect_allele.exposure,
#     A2 = other_allele.exposure,
#     freq = eaf.exposure,
#     b = beta.exposure,
#     se = se.exposure,
#     p = pval.exposure,
#     N = samplesize.exposure
#   )

# temp_input = tempfile()
# write_delim(dt, temp_input, num_threads = 1)

# cmd = glue("{gcta_path} --bfile {bfile} --cojo-file {temp_input} --cojo-slct --cojo-p 1e-3 --out {temp_input}")

# system(cmd)
