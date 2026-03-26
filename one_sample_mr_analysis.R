#!/usr/bin/env Rscript
###############################################################################
# One-Sample Mendelian Randomization (2SLS / 2SRI) Analysis
# Exposure: AIP / CHG / TyG (own GWAS, PC1-10, no BMI)
# Sample   : same individuals used for GWAS and individual-level MR
#
# This version removes dependency on MendelianRandomization package.
# Summary-level sensitivity (IVW / weighted median / MR-Egger / LOO / Q) is
# implemented with base R formulas.
###############################################################################

rm(list = ls()); gc()
options(stringsAsFactors = FALSE, scipen = 999)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required. Please install: install.packages('data.table')", call. = FALSE)
}
if (!requireNamespace("AER", quietly = TRUE)) {
  stop("Package 'AER' is required for ivreg. Please install: install.packages('AER')", call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(AER)   # ivreg
})

HAS_PRESSO <- requireNamespace("MRPRESSO", quietly = TRUE)
if (HAS_PRESSO) suppressPackageStartupMessages(library(MRPRESSO))

set.seed(20260326)

# =============================================================================
# 路径配置（按你的实际项目结构）
# =============================================================================
BASE        <- "/md3800f/file/Grade2024/mozhiwen/phd-T"
GWAS_DIR    <- file.path(BASE, "GWAS_PC10/gwas_primary")
CLUMP_DIR   <- file.path(BASE, "GWAS_PC10/clump")
PHENO_FILE  <- "/md3800f/file/Grade2024/mozhiwen/phd-T/GWAS2/3605人MT分析数据-表型加PCs版.txt"
BFILE       <- "/md3800f/file/Grade2024/mozhiwen/mt-staar/merge_6"
PCA_FILE    <- "/md3800f/file/Grade2024/mozhiwen/mt-staar/PCA_chrall_10PC.eigenvec"
OUT_DIR     <- file.path(BASE, "MR/OS_MR_PC10")
PLINK       <- Sys.which("plink")
if (PLINK == "") PLINK <- "plink"

for (d in c("results", "diagnostics", "intermediate", "instruments")) {
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# 需要分析的暴露和结局
# =============================================================================
exposures  <- c("AIP", "CHG", "TyG")
thresholds <- c("5e-8", "1e-6")

outcomes <- list(
  T2D = list(type = "binary",
             label = "Type 2 Diabetes",
             aliases = c("T2D", "t2d", "DM2", "dm2", "type2dm", "type2_diabetes", "diabetes")),
  HTN = list(type = "binary",
             label = "Hypertension",
             aliases = c("HTN", "htn", "hypertension", "high_bp", "highbp", "sbp_htn"))
)

# =============================================================================
# 工具函数
# =============================================================================
stop_if_missing <- function(x, msg) if (!x) stop(msg, call. = FALSE)

find_first_col <- function(nms, candidates, regex = FALSE) {
  if (!regex) {
    hit <- candidates[candidates %in% nms]
    return(if (length(hit)) hit[1] else NA_character_)
  }
  for (pat in candidates) {
    h <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(h)) return(h[1])
  }
  NA_character_
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

choose_best_id_transform <- function(pheno_id, geno_id) {
  geno_id <- as.character(geno_id)
  candidates <- list(
    identity = as.character(pheno_id),
    dup_self = paste0(pheno_id, "_", pheno_id),
    strip_dup = sub("^(\\d+)_\\1$", "\\1", as.character(pheno_id)),
    pad5 = suppressWarnings(sprintf("%05d", as.integer(pheno_id))),
    intchar = suppressWarnings(as.character(as.integer(pheno_id)))
  )
  overlaps <- sapply(candidates, function(v) sum(v %in% geno_id, na.rm = TRUE))
  best <- names(which.max(overlaps))[1]
  list(best = best, values = candidates[[best]], overlap = max(overlaps))
}

extract_coef <- function(fit, term) {
  cf <- summary(fit)$coefficients
  rn <- rownames(cf)
  idx <- which(rn == term)
  if (!length(idx)) return(list(beta = NA_real_, se = NA_real_, p = NA_real_, z = NA_real_))
  stat_col <- if (ncol(cf) >= 3) 3 else NA_integer_
  p_col    <- if (ncol(cf) >= 4) 4 else NA_integer_
  list(beta = unname(cf[idx, 1]),
       se   = unname(cf[idx, 2]),
       p    = if (!is.na(p_col)) unname(cf[idx, p_col]) else NA_real_,
       z    = if (!is.na(stat_col)) unname(cf[idx, stat_col]) else NA_real_)
}

first_stage_stats <- function(df, exposure_col, score_col, covars) {
  form0 <- as.formula(paste(exposure_col, "~", paste(covars, collapse = " + ")))
  form1 <- as.formula(paste(exposure_col, "~", paste(c(score_col, covars), collapse = " + ")))
  m0 <- lm(form0, data = df)
  m1 <- lm(form1, data = df)
  a  <- anova(m0, m1)
  rss0 <- deviance(m0)
  rss1 <- deviance(m1)
  partial_r2 <- (rss0 - rss1) / rss0
  data.table(
    score = score_col,
    n = nobs(m1),
    partial_F = a$F[2],
    p_first_stage = a$`Pr(>F)`[2],
    partial_R2 = partial_r2,
    beta_score = coef(summary(m1))[score_col, 1],
    se_score   = coef(summary(m1))[score_col, 2]
  )
}

run_2sri_binary <- function(df, outcome_col, exposure_col, score_col, covars) {
  form1 <- as.formula(paste(exposure_col, "~", paste(c(score_col, covars), collapse = " + ")))
  fs <- lm(form1, data = df)
  df$._xhat <- fitted(fs)
  df$._res1 <- residuals(fs)
  form2 <- as.formula(paste(outcome_col, "~ ._xhat + ._res1 +", paste(covars, collapse = " + ")))
  ss <- glm(form2, data = df, family = binomial())
  est <- extract_coef(ss, "._xhat")
  data.table(
    method = "2SRI_logistic",
    beta = est$beta,
    se   = est$se,
    pval = est$p,
    OR   = exp(est$beta),
    OR_lo = exp(est$beta - 1.96 * est$se),
    OR_hi = exp(est$beta + 1.96 * est$se)
  )
}

run_2sls <- function(df, outcome_col, exposure_col, score_col, covars, binary_outcome = FALSE) {
  rhs1 <- paste(c(exposure_col, covars), collapse = " + ")
  rhs2 <- paste(c(score_col, covars), collapse = " + ")
  form <- as.formula(paste(outcome_col, "~", rhs1, "|", rhs2))
  fit <- AER::ivreg(form, data = df)
  est <- extract_coef(fit, exposure_col)
  out <- data.table(
    method = if (binary_outcome) "2SLS_linear_probability" else "2SLS",
    beta = est$beta,
    se   = est$se,
    pval = est$p,
    OR = if (binary_outcome) NA_real_ else exp(est$beta),
    OR_lo = if (binary_outcome) NA_real_ else exp(est$beta - 1.96 * est$se),
    OR_hi = if (binary_outcome) NA_real_ else exp(est$beta + 1.96 * est$se)
  )
  diag <- tryCatch(summary(fit, diagnostics = TRUE)$diagnostics, error = function(e) NULL)
  list(result = out, diagnostics = diag)
}

run_observational <- function(df, outcome_col, exposure_col, covars, type) {
  form <- as.formula(paste(outcome_col, "~", paste(c(exposure_col, covars), collapse = " + ")))
  if (type == "binary") {
    fit <- glm(form, data = df, family = binomial())
    est <- extract_coef(fit, exposure_col)
    data.table(method = "Observational_logistic",
               beta = est$beta, se = est$se, pval = est$p,
               OR = exp(est$beta), OR_lo = exp(est$beta - 1.96 * est$se), OR_hi = exp(est$beta + 1.96 * est$se))
  } else {
    fit <- lm(form, data = df)
    est <- extract_coef(fit, exposure_col)
    data.table(method = "Observational_linear",
               beta = est$beta, se = est$se, pval = est$p,
               OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_)
  }
}

# -------- summary MR methods without MendelianRandomization --------
ivw_estimate <- function(beta_exp, se_exp, beta_out, se_out) {
  ratio <- beta_out / beta_exp
  se_ratio <- abs(se_out / beta_exp)
  w <- 1 / (se_ratio^2)
  b <- sum(w * ratio) / sum(w)
  se <- sqrt(1 / sum(w))
  p <- 2 * pnorm(-abs(b / se))
  list(beta = b, se = se, p = p, ratio = ratio, se_ratio = se_ratio, w = w)
}

weighted_median_estimate <- function(ratio, se_ratio) {
  w <- 1 / (se_ratio^2)
  ord <- order(ratio)
  r <- ratio[ord]; ww <- w[ord] / sum(w)
  cw <- cumsum(ww)
  idx <- which(cw >= 0.5)[1]
  b <- r[idx]
  # simple bootstrap SE
  set.seed(20260326)
  B <- 1000
  bs <- rep(NA_real_, B)
  n <- length(ratio)
  for (i in seq_len(B)) {
    id <- sample.int(n, n, replace = TRUE)
    rr <- ratio[id]; ss <- se_ratio[id]
    ww2 <- 1 / (ss^2)
    o2 <- order(rr)
    rr <- rr[o2]; ww2 <- ww2[o2] / sum(ww2)
    bs[i] <- rr[which(cumsum(ww2) >= 0.5)[1]]
  }
  se <- sd(bs, na.rm = TRUE)
  p <- ifelse(is.finite(se) && se > 0, 2 * pnorm(-abs(b / se)), NA_real_)
  list(beta = b, se = se, p = p)
}

egger_estimate <- function(beta_exp, se_exp, beta_out, se_out) {
  w <- 1 / (se_out^2)
  fit <- lm(beta_out ~ beta_exp, weights = w)
  s <- summary(fit)$coefficients
  if (!all(c("(Intercept)", "beta_exp") %in% rownames(s))) return(NULL)
  list(
    slope = s["beta_exp", 1],
    slope_se = s["beta_exp", 2],
    slope_p = s["beta_exp", 4],
    intercept = s["(Intercept)", 1],
    intercept_se = s["(Intercept)", 2],
    intercept_p = s["(Intercept)", 4]
  )
}

run_summary_mr <- function(dat, exposure, outcome, threshold, outcome_label) {
  dat <- dat[complete.cases(dat)]
  n_iv <- nrow(dat)
  stop_if_missing(n_iv >= 1, paste0("No valid SNP summary rows for ", exposure, " -> ", outcome))

  if (n_iv == 1) {
    beta <- dat$beta_out / dat$beta_exp
    se   <- dat$se_out / abs(dat$beta_exp)
    p    <- 2 * pnorm(-abs(beta / se))
    return(list(
      main = data.table(exposure = exposure, outcome = outcome, outcome_label = outcome_label,
                        threshold = threshold, n_iv = n_iv,
                        method = "Wald_ratio",
                        beta = beta, se = se, pval = p,
                        OR = exp(beta), OR_lo = exp(beta - 1.96 * se), OR_hi = exp(beta + 1.96 * se)),
      heterogeneity = NULL,
      egger_intercept = NULL,
      presso = NULL,
      loo = NULL
    ))
  }

  res <- list()
  ivw <- ivw_estimate(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out)
  res[["IVW"]] <- data.table(exposure = exposure, outcome = outcome, outcome_label = outcome_label,
                              threshold = threshold, n_iv = n_iv,
                              method = "IVW",
                              beta = ivw$beta, se = ivw$se, pval = ivw$p,
                              OR = exp(ivw$beta), OR_lo = exp(ivw$beta - 1.96 * ivw$se), OR_hi = exp(ivw$beta + 1.96 * ivw$se))

  Q <- sum(ivw$w * (ivw$ratio - ivw$beta)^2)
  heterogeneity <- data.table(exposure = exposure, outcome = outcome, outcome_label = outcome_label,
                              threshold = threshold, n_iv = n_iv,
                              Q = Q, Q_df = n_iv - 1, Q_pval = pchisq(Q, df = n_iv - 1, lower.tail = FALSE))

  if (n_iv >= 3) {
    wm <- weighted_median_estimate(ivw$ratio, ivw$se_ratio)
    res[["Weighted_median"]] <- data.table(exposure = exposure, outcome = outcome, outcome_label = outcome_label,
                                            threshold = threshold, n_iv = n_iv,
                                            method = "Weighted_median",
                                            beta = wm$beta, se = wm$se, pval = wm$p,
                                            OR = exp(wm$beta), OR_lo = exp(wm$beta - 1.96 * wm$se), OR_hi = exp(wm$beta + 1.96 * wm$se))
  }

  egger_intercept <- NULL
  if (n_iv >= 4) {
    eg <- egger_estimate(dat$beta_exp, dat$se_exp, dat$beta_out, dat$se_out)
    if (!is.null(eg)) {
      res[["MR_Egger"]] <- data.table(exposure = exposure, outcome = outcome, outcome_label = outcome_label,
                                       threshold = threshold, n_iv = n_iv,
                                       method = "MR_Egger",
                                       beta = eg$slope, se = eg$slope_se, pval = eg$slope_p,
                                       OR = exp(eg$slope), OR_lo = exp(eg$slope - 1.96 * eg$slope_se), OR_hi = exp(eg$slope + 1.96 * eg$slope_se))
      egger_intercept <- data.table(exposure = exposure, outcome = outcome, outcome_label = outcome_label,
                                    threshold = threshold, n_iv = n_iv,
                                    intercept = eg$intercept, se = eg$intercept_se, pval = eg$intercept_p)
    }
  }

  presso_res <- NULL
  if (n_iv >= 4 && HAS_PRESSO) {
    pr <- tryCatch(
      MRPRESSO::mr_presso(BetaOutcome = "beta_out", BetaExposure = "beta_exp",
                          SdOutcome = "se_out", SdExposure = "se_exp",
                          data = as.data.frame(dat),
                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                          NbDistribution = 3000, SignifThreshold = 0.05),
      error = function(e) NULL
    )
    if (!is.null(pr)) {
      gp <- tryCatch(pr$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA_real_)
      ce <- tryCatch(pr$`Main MR results`[1, "Causal Estimate"], error = function(e) NA_real_)
      se <- tryCatch(pr$`Main MR results`[1, "Sd"], error = function(e) NA_real_)
      pv <- tryCatch(pr$`Main MR results`[1, "P-value"], error = function(e) NA_real_)
      presso_res <- data.table(exposure = exposure, outcome = outcome, outcome_label = outcome_label,
                               threshold = threshold, n_iv = n_iv,
                               method = "MR_PRESSO",
                               beta = ce, se = se, pval = pv,
                               OR = exp(ce), OR_lo = exp(ce - 1.96 * se), OR_hi = exp(ce + 1.96 * se),
                               presso_global_p = gp)
    }
  }

  loo <- NULL
  if (n_iv >= 3) {
    tmp <- vector("list", n_iv)
    for (i in seq_len(n_iv)) {
      sub <- dat[-i]
      if (nrow(sub) == 1) {
        b <- sub$beta_out / sub$beta_exp
        s <- sub$se_out / abs(sub$beta_exp)
        p <- 2 * pnorm(-abs(b / s))
      } else {
        o <- ivw_estimate(sub$beta_exp, sub$se_exp, sub$beta_out, sub$se_out)
        b <- o$beta; s <- o$se; p <- o$p
      }
      tmp[[i]] <- data.table(excluded_snp = dat$SNP[i], beta = b, se = s, pval = p,
                             OR = exp(b), OR_lo = exp(b - 1.96 * s), OR_hi = exp(b + 1.96 * s))
    }
    loo <- rbindlist(tmp)
    loo[, `:=`(exposure = exposure, outcome = outcome, threshold = threshold)]
  }

  list(main = rbindlist(res, fill = TRUE),
       heterogeneity = heterogeneity,
       egger_intercept = egger_intercept,
       presso = presso_res,
       loo = loo)
}

build_crossfit_score <- function(df, snp_cols, exposure_col, covars, K = 5) {
  n <- nrow(df)
  fold_id <- sample(rep(seq_len(K), length.out = n))
  cf_score <- rep(NA_real_, n)
  betamat <- matrix(NA_real_, nrow = K, ncol = length(snp_cols), dimnames = list(paste0("Fold", 1:K), snp_cols))

  for (k in seq_len(K)) {
    idx_te <- which(fold_id == k)
    idx_tr <- setdiff(seq_len(n), idx_te)
    dtr <- df[idx_tr]
    dte <- df[idx_te]

    b <- rep(NA_real_, length(snp_cols))
    names(b) <- snp_cols
    for (s in snp_cols) {
      fm <- lm(as.formula(paste(exposure_col, "~", paste(c(s, covars), collapse = " + "))), data = dtr)
      est <- extract_coef(fm, s)
      b[s] <- est$beta
    }
    betamat[k, ] <- b
    cf_score[idx_te] <- as.numeric(as.matrix(dte[, ..snp_cols]) %*% b)
  }

  list(score = cf_score,
       fold = fold_id,
       beta_matrix = as.data.table(betamat, keep.rownames = "fold"))
}

read_pheno <- function() {
  stop_if_missing(file.exists(PHENO_FILE), paste0("Phenotype file not found: ", PHENO_FILE))
  ph <- data.table::fread(PHENO_FILE)

  iid_col <- find_first_col(names(ph), c("IID", "iid", "sample.id", "sample_id", "ID", "id"))
  stop_if_missing(!is.na(iid_col), "Cannot find IID/sample ID column in phenotype file.")
  setnames(ph, iid_col, "IID_PHENO")
  ph[, IID_PHENO := as.character(IID_PHENO)]

  age_col <- find_first_col(names(ph), c("age", "Age", "AGE", "age_years"))
  sex_col <- find_first_col(names(ph), c("sex", "Sex", "SEX", "gender", "Gender", "GENDER"))
  stop_if_missing(!is.na(age_col), "Cannot find age column in phenotype file.")
  stop_if_missing(!is.na(sex_col), "Cannot find sex/gender column in phenotype file.")

  if (age_col != "age") setnames(ph, age_col, "age")
  if (sex_col != "sex") setnames(ph, sex_col, "sex")

  for (e in exposures) stop_if_missing(e %in% names(ph), paste0("Exposure column missing in phenotype file: ", e))

  outcome_map <- list()
  for (nm in names(outcomes)) {
    hit <- find_first_col(names(ph), outcomes[[nm]]$aliases)
    stop_if_missing(!is.na(hit), paste0("Cannot find outcome column for ", nm,
                                       ". Checked aliases: ", paste(outcomes[[nm]]$aliases, collapse = ", ")))
    outcome_map[[nm]] <- hit
  }

  need_pc_merge <- !all(paste0("PC", 1:10) %in% names(ph))
  if (need_pc_merge) {
    stop_if_missing(file.exists(PCA_FILE), paste0("PC file not found: ", PCA_FILE))
    pca <- data.table::fread(PCA_FILE, header = FALSE)
    stop_if_missing(ncol(pca) >= 12, "PCA file must contain FID IID PC1-PC10.")
    setnames(pca, 1:12, c("FID", "IID", paste0("PC", 1:10)))
    pca[, IID := as.character(IID)]

    best <- choose_best_id_transform(ph$IID_PHENO, pca$IID)
    ph[, IID := best$values]
    stop_if_missing(best$overlap > 0, "Phenotype IDs cannot be matched to PCA IDs.")
    ph <- merge(ph, pca[, c("IID", paste0("PC", 1:10)), with = FALSE], by = "IID", all.x = TRUE)
  }

  if (!"IID" %in% names(ph)) ph[, IID := as.character(IID_PHENO)]

  ph[, sex := as.character(sex)]
  ph[sex %in% c("M", "Male", "male", "1"), sex := "1"]
  ph[sex %in% c("F", "Female", "female", "2"), sex := "2"]
  ph[, sex := safe_numeric(sex)]
  ph[, age := safe_numeric(age)]

  for (nm in names(outcome_map)) {
    col <- outcome_map[[nm]]
    ph[[col]] <- safe_numeric(ph[[col]])
  }

  for (e in exposures) ph[[e]] <- safe_numeric(ph[[e]])
  for (pc in paste0("PC", 1:10)) ph[[pc]] <- safe_numeric(ph[[pc]])

  list(data = ph, outcome_map = outcome_map)
}

read_iv <- function(exposure, threshold) {
  snp_file <- if (threshold == "5e-8") file.path(CLUMP_DIR, sprintf("%s_lead_snps.txt", exposure))
  else file.path(CLUMP_DIR, sprintf("%s_lead_snps_1e6.txt", exposure))
  stop_if_missing(file.exists(snp_file), paste0("IV file not found: ", snp_file))
  lead <- trimws(readLines(snp_file))
  lead <- lead[lead != ""]
  stop_if_missing(length(lead) > 0, paste0("No lead SNPs in ", snp_file))

  gwas_file <- file.path(GWAS_DIR, sprintf("%s.%s.glm.linear", exposure, exposure))
  stop_if_missing(file.exists(gwas_file), paste0("GWAS file not found: ", gwas_file))
  gwas <- data.table::fread(gwas_file)
  if ("#CHROM" %in% names(gwas)) setnames(gwas, "#CHROM", "CHR")

  iv <- gwas[ID %in% lead, .(
    SNP = ID,
    CHR = CHR,
    POS = POS,
    effect_allele = A1,
    other_allele  = ifelse(A1 == REF, ALT, REF),
    beta_gwas = BETA,
    se_gwas   = SE,
    pval_gwas = P,
    eaf = A1_FREQ
  )]
  stop_if_missing(nrow(iv) > 0, paste0("No IVs matched in GWAS for ", exposure, " @ ", threshold))

  iv[, F_gwas := (beta_gwas / se_gwas)^2]
  iv[, R2_gwas := (beta_gwas^2) / (beta_gwas^2 + se_gwas^2 * 3605)]
  iv[]
}

extract_genotypes <- function(iv, exposure, threshold) {
  snp_file <- file.path(OUT_DIR, "intermediate", sprintf("%s_%s_snps.txt", exposure, threshold))
  fwrite(iv[, .(SNP)], snp_file, col.names = FALSE)

  outprefix <- file.path(OUT_DIR, "intermediate", sprintf("%s_%s_extract", exposure, threshold))
  raw_file  <- paste0(outprefix, ".raw")
  if (!file.exists(raw_file)) {
    cmd <- sprintf("%s --bfile %s --extract %s --recode A --out %s",
                   shQuote(PLINK), shQuote(BFILE), shQuote(snp_file), shQuote(outprefix))
    cat("Running:", cmd, "\n")
    rc <- system(cmd)
    stop_if_missing(rc == 0 && file.exists(raw_file), paste0("PLINK extraction failed for ", exposure, " @ ", threshold))
  }

  gt <- data.table::fread(raw_file)
  gt[, FID := as.character(FID)]
  gt[, IID := as.character(IID)]

  meta_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  gcols <- setdiff(names(gt), meta_cols)
  stop_if_missing(length(gcols) > 0, paste0("No genotype columns found in ", raw_file))

  parse_base <- function(x) sub("_([^_]*)$", "", x)
  parse_allele <- function(x) sub("^.*_", "", x)

  mapped <- list()
  for (i in seq_len(nrow(iv))) {
    snp <- iv$SNP[i]
    ea  <- toupper(iv$effect_allele[i])
    oa  <- toupper(iv$other_allele[i])
    hit <- gcols[parse_base(gcols) == snp]
    stop_if_missing(length(hit) == 1, paste0("Cannot uniquely find genotype column for SNP: ", snp, " in ", raw_file))
    col <- hit[1]
    counted <- toupper(parse_allele(col))
    dose <- gt[[col]]
    if (counted == ea) dose_ea <- dose else if (counted == oa) dose_ea <- 2 - dose else {
      stop(paste0("Counted allele mismatch for ", snp, ": counted=", counted,
                  ", effect=", ea, ", other=", oa), call. = FALSE)
    }
    mapped[[snp]] <- dose_ea
  }

  geno <- data.table(IID = gt$IID)
  for (s in names(mapped)) geno[[s]] <- mapped[[s]]
  geno
}

merge_pheno_geno <- function(ph, geno) {
  best <- choose_best_id_transform(ph$IID, geno$IID)
  stop_if_missing(best$overlap > 0, "Phenotype IDs cannot be matched to genotype IDs.")
  ph2 <- copy(ph)
  ph2[, IID_M := best$values]
  mg <- merge(ph2, geno, by.x = "IID_M", by.y = "IID", all = FALSE)
  mg[, IID := IID_M]
  mg[, IID_M := NULL]
  mg
}

make_same_sample_summary <- function(df, snps, exposure_col, outcome_col, outcome_type, covars) {
  out <- vector("list", length(snps))
  for (i in seq_along(snps)) {
    s <- snps[i]
    fm_exp <- lm(as.formula(paste(exposure_col, "~", paste(c(s, covars), collapse = " + "))), data = df)
    e <- extract_coef(fm_exp, s)

    if (outcome_type == "binary") {
      fm_out <- glm(as.formula(paste(outcome_col, "~", paste(c(s, covars), collapse = " + "))), data = df, family = binomial())
    } else {
      fm_out <- lm(as.formula(paste(outcome_col, "~", paste(c(s, covars), collapse = " + "))), data = df)
    }
    o <- extract_coef(fm_out, s)

    out[[i]] <- data.table(
      SNP = s,
      beta_exp = e$beta, se_exp = e$se, pval_exp = e$p,
      F_exp = ifelse(is.na(e$se) || e$se == 0, NA_real_, (e$beta / e$se)^2),
      beta_out = o$beta, se_out = o$se, pval_out = o$p
    )
  }
  rbindlist(out)
}

run_one_combo <- function(ph, outcome_map, exposure, threshold, outcome_nm) {
  cat(sprintf("\n===== %s -> %s @ %s =====\n", exposure, outcome_nm, threshold))

  iv <- read_iv(exposure, threshold)
  fwrite(iv, file.path(OUT_DIR, "instruments", sprintf("%s_%s_IVs_from_GWAS.csv", exposure, threshold)))

  geno <- extract_genotypes(iv, exposure, threshold)
  dat  <- merge_pheno_geno(ph, geno)

  outcome_col  <- outcome_map[[outcome_nm]]
  outcome_type <- outcomes[[outcome_nm]]$type
  outcome_lab  <- outcomes[[outcome_nm]]$label
  covars <- c("age", "sex", paste0("PC", 1:10))
  req <- c(exposure, outcome_col, covars, iv$SNP)
  dat <- dat[complete.cases(dat[, ..req])]
  stop_if_missing(nrow(dat) > 100, paste0("Too few complete cases for ", exposure, " -> ", outcome_nm))

  snps <- iv$SNP
  w <- iv$beta_gwas; names(w) <- iv$SNP
  dat[, GRS_weighted := as.numeric(as.matrix(dat[, ..snps]) %*% w)]
  dat[, GRS_unweighted := rowSums(.SD), .SDcols = snps]

  obs <- run_observational(dat, outcome_col, exposure, covars, outcome_type)

  fs_main <- first_stage_stats(dat, exposure, "GRS_weighted", covars)
  fs_unw  <- first_stage_stats(dat, exposure, "GRS_unweighted", covars)
  fs_main[, score_type := "weighted_GRS"]
  fs_unw[,  score_type := "unweighted_GRS"]
  fs_all <- rbind(fs_main, fs_unw, fill = TRUE)

  if (outcome_type == "binary") {
    main_res <- run_2sri_binary(dat, outcome_col, exposure, "GRS_weighted", covars)
    lpm <- run_2sls(dat, outcome_col, exposure, "GRS_weighted", covars, binary_outcome = TRUE)
    res_main <- rbindlist(list(
      cbind(data.table(exposure = exposure, outcome = outcome_nm, outcome_label = outcome_lab,
                       threshold = threshold, n = nrow(dat), n_iv = length(snps), score = "weighted_GRS"), main_res),
      cbind(data.table(exposure = exposure, outcome = outcome_nm, outcome_label = outcome_lab,
                       threshold = threshold, n = nrow(dat), n_iv = length(snps), score = "weighted_GRS"), lpm$result),
      cbind(data.table(exposure = exposure, outcome = outcome_nm, outcome_label = outcome_lab,
                       threshold = threshold, n = nrow(dat), n_iv = length(snps), score = "NA"), obs)
    ), fill = TRUE)
    diag_ivreg <- lpm$diagnostics
  } else {
    tsls <- run_2sls(dat, outcome_col, exposure, "GRS_weighted", covars, binary_outcome = FALSE)
    res_main <- rbindlist(list(
      cbind(data.table(exposure = exposure, outcome = outcome_nm, outcome_label = outcome_lab,
                       threshold = threshold, n = nrow(dat), n_iv = length(snps), score = "weighted_GRS"), tsls$result),
      cbind(data.table(exposure = exposure, outcome = outcome_nm, outcome_label = outcome_lab,
                       threshold = threshold, n = nrow(dat), n_iv = length(snps), score = "NA"), obs)
    ), fill = TRUE)
    diag_ivreg <- tsls$diagnostics
  }

  smry <- make_same_sample_summary(dat, snps, exposure, outcome_col, outcome_type, covars)
  smry <- merge(iv[, .(SNP, effect_allele, other_allele, beta_gwas, se_gwas, F_gwas, R2_gwas)], smry, by = "SNP", all.x = TRUE)
  fwrite(smry, file.path(OUT_DIR, "diagnostics", sprintf("%s_%s_%s_same_sample_SNP_associations.csv", exposure, outcome_nm, threshold)))

  mrdiag <- run_summary_mr(smry[, .(SNP, beta_exp, se_exp, beta_out, se_out)], exposure, outcome_nm, threshold, outcome_lab)

  cf <- build_crossfit_score(dat, snps, exposure, covars, K = 5)
  dat[, GRS_crossfit := cf$score]
  fs_cf <- first_stage_stats(dat, exposure, "GRS_crossfit", covars)
  fs_cf[, score_type := "crossfit_weighted_GRS"]
  fs_all <- rbind(fs_all, fs_cf, fill = TRUE)

  if (outcome_type == "binary") {
    cf_res <- run_2sri_binary(dat, outcome_col, exposure, "GRS_crossfit", covars)
    cf_lpm <- run_2sls(dat, outcome_col, exposure, "GRS_crossfit", covars, binary_outcome = TRUE)
    res_cf <- rbindlist(list(
      cbind(data.table(exposure = exposure, outcome = outcome_nm, outcome_label = outcome_lab,
                       threshold = threshold, n = nrow(dat), n_iv = length(snps), score = "crossfit_weighted_GRS"), cf_res),
      cbind(data.table(exposure = exposure, outcome = outcome_nm, outcome_label = outcome_lab,
                       threshold = threshold, n = nrow(dat), n_iv = length(snps), score = "crossfit_weighted_GRS"), cf_lpm$result)
    ), fill = TRUE)
    diag_cf <- cf_lpm$diagnostics
  } else {
    cf_tsls <- run_2sls(dat, outcome_col, exposure, "GRS_crossfit", covars, binary_outcome = FALSE)
    res_cf <- cbind(data.table(exposure = exposure, outcome = outcome_nm, outcome_label = outcome_lab,
                               threshold = threshold, n = nrow(dat), n_iv = length(snps), score = "crossfit_weighted_GRS"), cf_tsls$result)
    diag_cf <- cf_tsls$diagnostics
  }

  fwrite(fs_all, file.path(OUT_DIR, "diagnostics", sprintf("%s_%s_%s_first_stage.csv", exposure, outcome_nm, threshold)))
  fwrite(cf$beta_matrix, file.path(OUT_DIR, "diagnostics", sprintf("%s_%s_%s_crossfit_training_betas.csv", exposure, outcome_nm, threshold)))

  if (!is.null(mrdiag$loo)) fwrite(mrdiag$loo, file.path(OUT_DIR, "diagnostics", sprintf("%s_%s_%s_leave_one_out.csv", exposure, outcome_nm, threshold)))
  if (!is.null(mrdiag$heterogeneity)) fwrite(mrdiag$heterogeneity, file.path(OUT_DIR, "diagnostics", sprintf("%s_%s_%s_heterogeneity.csv", exposure, outcome_nm, threshold)))
  if (!is.null(mrdiag$egger_intercept)) fwrite(mrdiag$egger_intercept, file.path(OUT_DIR, "diagnostics", sprintf("%s_%s_%s_egger_intercept.csv", exposure, outcome_nm, threshold)))
  if (!is.null(mrdiag$presso)) fwrite(mrdiag$presso, file.path(OUT_DIR, "diagnostics", sprintf("%s_%s_%s_mrpresso.csv", exposure, outcome_nm, threshold)))

  list(
    main_results = rbindlist(list(res_main, res_cf), fill = TRUE),
    first_stage = fs_all,
    snp_assoc = smry,
    mrdiag_main = mrdiag$main,
    mrdiag_heterogeneity = mrdiag$heterogeneity,
    mrdiag_egger = mrdiag$egger_intercept,
    mrdiag_presso = mrdiag$presso,
    ivreg_diag_main = if (!is.null(diag_ivreg)) as.data.table(diag_ivreg, keep.rownames = "test") else NULL,
    ivreg_diag_crossfit = if (!is.null(diag_cf)) as.data.table(diag_cf, keep.rownames = "test") else NULL
  )
}

cat("###################################################################\n")
cat("# One-Sample MR / 2SLS / 2SRI\n")
cat("#", format(Sys.time()), "\n")
cat("###################################################################\n\n")
cat("BASE       :", BASE, "\n")
cat("PHENO_FILE :", PHENO_FILE, "\n")
cat("BFILE      :", BFILE, "\n")
cat("GWAS_DIR   :", GWAS_DIR, "\n")
cat("CLUMP_DIR  :", CLUMP_DIR, "\n")
cat("OUT_DIR    :", OUT_DIR, "\n\n")
if (!HAS_PRESSO) cat("NOTE: MRPRESSO not installed; MR-PRESSO step will be skipped.\n\n")

obj <- read_pheno()
ph  <- obj$data
outcome_map <- obj$outcome_map

all_main <- list(); all_first <- list(); all_diag_main <- list(); all_het <- list()
all_egger <- list(); all_presso <- list(); all_ivreg1 <- list(); all_ivreg2 <- list()

idx <- 0
for (e in exposures) {
  for (th in thresholds) {
    for (o in names(outcomes)) {
      idx <- idx + 1
      ans <- run_one_combo(ph, outcome_map, e, th, o)
      all_main[[idx]]  <- ans$main_results
      all_first[[idx]] <- ans$first_stage
      all_diag_main[[idx]] <- ans$mrdiag_main
      all_het[[idx]]   <- ans$mrdiag_heterogeneity
      all_egger[[idx]] <- ans$mrdiag_egger
      all_presso[[idx]] <- ans$mrdiag_presso
      all_ivreg1[[idx]] <- ans$ivreg_diag_main
      all_ivreg2[[idx]] <- ans$ivreg_diag_crossfit
    }
  }
}

main_tab  <- rbindlist(all_main, fill = TRUE)
first_tab <- rbindlist(all_first, fill = TRUE)
mr_tab    <- rbindlist(all_diag_main, fill = TRUE)
bind_or_empty <- function(x) {
  y <- Filter(Negate(is.null), x)
  if (!length(y)) return(data.table())
  rbindlist(y, fill = TRUE)
}

het_tab <- bind_or_empty(all_het); eg_tab <- bind_or_empty(all_egger); pr_tab <- bind_or_empty(all_presso)
iv1_tab <- bind_or_empty(all_ivreg1); iv2_tab <- bind_or_empty(all_ivreg2)

fwrite(main_tab,  file.path(OUT_DIR, "results", "one_sample_MR_main_results.csv"))
fwrite(first_tab, file.path(OUT_DIR, "results", "one_sample_MR_first_stage_summary.csv"))
fwrite(mr_tab,    file.path(OUT_DIR, "results", "one_sample_MR_summary_based_sensitivity.csv"))
if (nrow(het_tab)) fwrite(het_tab, file.path(OUT_DIR, "results", "one_sample_MR_heterogeneity.csv"))
if (nrow(eg_tab))  fwrite(eg_tab,  file.path(OUT_DIR, "results", "one_sample_MR_egger_intercept.csv"))
if (nrow(pr_tab))  fwrite(pr_tab,  file.path(OUT_DIR, "results", "one_sample_MR_mrpresso.csv"))
if (nrow(iv1_tab)) fwrite(iv1_tab, file.path(OUT_DIR, "results", "one_sample_MR_ivreg_diagnostics_main.csv"))
if (nrow(iv2_tab)) fwrite(iv2_tab, file.path(OUT_DIR, "results", "one_sample_MR_ivreg_diagnostics_crossfit.csv"))

cat("\n================== Summary ==================\n")
print(main_tab)
cat("\nOutput directory:\n", OUT_DIR, "\n")
cat("Done at:", format(Sys.time()), "\n")
