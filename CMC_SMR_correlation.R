## ============================================================
## CMCinc vs beta modulation analysis
##
## Methods:
## 1. Spearman correlations
## 2. Nonparametric bootstrap 95% CI (10,000 resamples)
## 3. Two-sided permutation test (10,000 permutations)
## 4. Partial Spearman correlations:
##      |MRBDmax| vs CMCinc controlling SRBSmax
##      SRBSmax vs CMCinc controlling |MRBDmax|
## 5. Theil-Sen robust trend lines for visualization
##
## Variables in original data:
## CMC = CMCinc
## MRBD = MRBDmax (negative)
## SRBS = SRBSmax (positive)
## ============================================================


## ------------------------------------------------------------
## 0) Load data
## ------------------------------------------------------------

df <- read.csv(
  "/Users/Desktop/correlation.csv",
  stringsAsFactors = FALSE
)


## ------------------------------------------------------------
## 1) Packages
## ------------------------------------------------------------

pkgs <- c(
  "ggplot2",
  "mblm",
  "dplyr"
)

to_install <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(to_install) > 0) {
  install.packages(to_install)
}

library(ggplot2)
library(mblm)
library(dplyr)


## ------------------------------------------------------------
## 2) Prepare variables
## ------------------------------------------------------------

stopifnot(
  all(c("CMC", "MRBD", "SRBS") %in% names(df))
)

df0 <- df %>%
  transmute(
    
    ## CMC increase
    CMCinc = as.numeric(CMC),
    
    ## MRBD is negative; use absolute magnitude
    MRBD_mag = abs(as.numeric(MRBD)),
    
    ## SRBS is positive
    SRBSmax = as.numeric(SRBS)
    
  ) %>%
  
  filter(
    is.finite(CMCinc),
    is.finite(MRBD_mag),
    is.finite(SRBSmax)
  ) %>%
  
  mutate(
    
    ## ΔSRBS-MRBD
    ## = SRBSmax - MRBDmax
    ## = SRBSmax + |MRBDmax| because MRBDmax is negative
    Delta_SRBS_MRBD = SRBSmax + MRBD_mag
  )


n <- nrow(df0)

if (n < 8) {
  stop("Too few complete cases.")
}

cat("\n============================================\n")
cat("Sample size\n")
cat("============================================\n")
cat("N =", n, "\n")


## ============================================================
## 3) Spearman correlation
##    + bootstrap 95% CI
##    + two-sided permutation test
## ============================================================

spearman_full <- function(
    x,
    y,
    R_boot = 10000,
    R_perm = 10000
) {
  
  ## complete finite cases
  ok <- is.finite(x) & is.finite(y)
  
  x <- x[ok]
  y <- y[ok]
  
  n <- length(x)
  
  if (n < 3) {
    stop("Too few observations.")
  }
  
  if (
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    stop("One variable has no variation.")
  }
  
  
  ## ----------------------------------------------------------
  ## Observed Spearman rho
  ## ----------------------------------------------------------
  
  rho_obs <- suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )
  
  
  ## ----------------------------------------------------------
  ## Nonparametric bootstrap 95% CI
  ## Resample subject pairs with replacement
  ## ----------------------------------------------------------
  
  rho_boot <- numeric(R_boot)
  
  for (i in seq_len(R_boot)) {
    
    idx <- sample.int(
      n,
      size = n,
      replace = TRUE
    )
    
    rho_boot[i] <- suppressWarnings(
      cor(
        x[idx],
        y[idx],
        method = "spearman"
      )
    )
  }
  
  
  CI <- quantile(
    rho_boot,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  
  ## ----------------------------------------------------------
  ## Two-sided permutation test
  ## ----------------------------------------------------------
  
  rho_perm <- numeric(R_perm)
  
  for (i in seq_len(R_perm)) {
    
    y_perm <- sample(
      y,
      size = n,
      replace = FALSE
    )
    
    rho_perm[i] <- suppressWarnings(
      cor(
        x,
        y_perm,
        method = "spearman"
      )
    )
  }
  
  
  ## finite-sample corrected permutation P
  p_perm <- (
    sum(
      abs(rho_perm) >= abs(rho_obs),
      na.rm = TRUE
    ) + 1
  ) / (R_perm + 1)
  
  
  ## ----------------------------------------------------------
  ## Return
  ## ----------------------------------------------------------
  
  list(
    rho = rho_obs,
    CI_low = CI[1],
    CI_high = CI[2],
    p_perm = p_perm,
    N = n,
    bootstrap_distribution = rho_boot,
    permutation_distribution = rho_perm
  )
}


## ============================================================
## 4) Partial Spearman correlation
##    + two-sided permutation test
##
## No bootstrap CI, matching the described Methods.
##
## Freedman-Lane residual permutation is used so that
## the effect of the control variable is retained.
## ============================================================

partial_spearman_perm <- function(
    x,
    y,
    z,
    R_perm = 10000
) {
  
  ## complete finite cases
  ok <- (
    is.finite(x) &
      is.finite(y) &
      is.finite(z)
  )
  
  x <- x[ok]
  y <- y[ok]
  z <- z[ok]
  
  n <- length(x)
  
  
  if (n < 4) {
    stop("Too few observations.")
  }
  
  if (
    length(unique(x)) < 2 ||
    length(unique(y)) < 2 ||
    length(unique(z)) < 2
  ) {
    stop("One variable has no variation.")
  }
  
  
  ## ----------------------------------------------------------
  ## Rank-transform variables
  ## Partial Spearman = partial Pearson correlation of ranks
  ## ----------------------------------------------------------
  
  rx <- rank(
    x,
    ties.method = "average"
  )
  
  ry <- rank(
    y,
    ties.method = "average"
  )
  
  rz <- rank(
    z,
    ties.method = "average"
  )
  
  
  ## ----------------------------------------------------------
  ## Residualize X and Y with respect to Z
  ## ----------------------------------------------------------
  
  fit_x <- lm(rx ~ rz)
  fit_y <- lm(ry ~ rz)
  
  x_res <- residuals(fit_x)
  y_res <- residuals(fit_y)
  
  
  ## observed partial Spearman rho
  rho_obs <- cor(
    x_res,
    y_res,
    method = "pearson"
  )
  
  
  ## ----------------------------------------------------------
  ## Freedman-Lane permutation
  ##
  ## Reduced model:
  ## ranked Y ~ ranked Z
  ##
  ## Permute residuals from this model while retaining
  ## the fitted contribution of Z.
  ## ----------------------------------------------------------
  
  y_fitted <- fitted(fit_y)
  y_residuals <- residuals(fit_y)
  
  rho_perm <- numeric(R_perm)
  
  
  for (i in seq_len(R_perm)) {
    
    ## permute reduced-model residuals
    perm_res <- sample(
      y_residuals,
      size = n,
      replace = FALSE
    )
    
    ## construct permuted ranked outcome
    ry_perm <- y_fitted + perm_res
    
    ## remove Z again
    y_perm_res <- residuals(
      lm(ry_perm ~ rz)
    )
    
    ## partial correlation under permutation
    rho_perm[i] <- cor(
      x_res,
      y_perm_res,
      method = "pearson"
    )
  }
  
  
  ## ----------------------------------------------------------
  ## Two-sided permutation P
  ## ----------------------------------------------------------
  
  p_perm <- (
    sum(
      abs(rho_perm) >= abs(rho_obs),
      na.rm = TRUE
    ) + 1
  ) / (R_perm + 1)
  
  
  ## ----------------------------------------------------------
  ## Return
  ## ----------------------------------------------------------
  
  list(
    rho = rho_obs,
    p_perm = p_perm,
    N = n,
    permutation_distribution = rho_perm
  )
}


## ============================================================
## 5) Run analyses
## ============================================================

set.seed(2026)

R_BOOT <- 10000
R_PERM <- 10000


## ------------------------------------------------------------
## 5.1 Ordinary Spearman correlations
## ------------------------------------------------------------

## |MRBDmax| vs CMCinc
sp_MRBD <- spearman_full(
  x = df0$MRBD_mag,
  y = df0$CMCinc,
  R_boot = R_BOOT,
  R_perm = R_PERM
)


## SRBSmax vs CMCinc
sp_SRBS <- spearman_full(
  x = df0$SRBSmax,
  y = df0$CMCinc,
  R_boot = R_BOOT,
  R_perm = R_PERM
)


## ΔSRBS-MRBD vs CMCinc
sp_DELTA <- spearman_full(
  x = df0$Delta_SRBS_MRBD,
  y = df0$CMCinc,
  R_boot = R_BOOT,
  R_perm = R_PERM
)


## ------------------------------------------------------------
## 5.2 Partial Spearman correlations
## ------------------------------------------------------------

## |MRBDmax| vs CMCinc
## controlling for SRBSmax
ps_MRBD <- partial_spearman_perm(
  x = df0$MRBD_mag,
  y = df0$CMCinc,
  z = df0$SRBSmax,
  R_perm = R_PERM
)


## SRBSmax vs CMCinc
## controlling for |MRBDmax|
ps_SRBS <- partial_spearman_perm(
  x = df0$SRBSmax,
  y = df0$CMCinc,
  z = df0$MRBD_mag,
  R_perm = R_PERM
)


## ============================================================
## 6) Ordinary Spearman results table
##
## These results contain:
## rho
## bootstrap 95% CI
## permutation P
## ============================================================

results_spearman <- data.frame(
  
  Analysis = c(
    "CMCinc vs |MRBDmax|",
    "CMCinc vs SRBSmax",
    "CMCinc vs DeltaSRBS-MRBD"
  ),
  
  rho = c(
    sp_MRBD$rho,
    sp_SRBS$rho,
    sp_DELTA$rho
  ),
  
  CI_2.5 = c(
    sp_MRBD$CI_low,
    sp_SRBS$CI_low,
    sp_DELTA$CI_low
  ),
  
  CI_97.5 = c(
    sp_MRBD$CI_high,
    sp_SRBS$CI_high,
    sp_DELTA$CI_high
  ),
  
  permutation_p = c(
    sp_MRBD$p_perm,
    sp_SRBS$p_perm,
    sp_DELTA$p_perm
  ),
  
  N = c(
    sp_MRBD$N,
    sp_SRBS$N,
    sp_DELTA$N
  )
)


cat("\n============================================\n")
cat("SPEARMAN CORRELATIONS\n")
cat("95% CI = 10,000-resample nonparametric bootstrap\n")
cat("P = two-sided 10,000-permutation test\n")
cat("============================================\n\n")

print(
  results_spearman,
  digits = 5,
  row.names = FALSE
)


## ============================================================
## 7) Partial Spearman results table
##
## No bootstrap CI
## ============================================================

results_partial <- data.frame(
  
  Analysis = c(
    "CMCinc vs |MRBDmax| controlling SRBSmax",
    "CMCinc vs SRBSmax controlling |MRBDmax|"
  ),
  
  partial_rho = c(
    ps_MRBD$rho,
    ps_SRBS$rho
  ),
  
  permutation_p = c(
    ps_MRBD$p_perm,
    ps_SRBS$p_perm
  ),
  
  N = c(
    ps_MRBD$N,
    ps_SRBS$N
  )
)


cat("\n============================================\n")
cat("PARTIAL SPEARMAN CORRELATIONS\n")
cat("P = two-sided 10,000-permutation test\n")
cat("============================================\n\n")

print(
  results_partial,
  digits = 5,
  row.names = FALSE
)


## ============================================================
## 8) Theil-Sen robust trend lines
##    Visualization only
## ============================================================

## mblm works most reliably with a base data.frame

df1 <- as.data.frame(df0)

df1$CMCinc <- as.numeric(
  df1$CMCinc
)

df1$MRBD_mag <- as.numeric(
  df1$MRBD_mag
)

df1$SRBSmax <- as.numeric(
  df1$SRBSmax
)

df1$Delta_SRBS_MRBD <- as.numeric(
  df1$Delta_SRBS_MRBD
)


## ------------------------------------------------------------
## Theil-Sen models
## ------------------------------------------------------------

ts_MRBD <- mblm::mblm(
  CMCinc ~ MRBD_mag,
  data = df1,
  repeated = FALSE
)


ts_SRBS <- mblm::mblm(
  CMCinc ~ SRBSmax,
  data = df1,
  repeated = FALSE
)


ts_DELTA <- mblm::mblm(
  CMCinc ~ Delta_SRBS_MRBD,
  data = df1,
  repeated = FALSE
)


cat("\n============================================\n")
cat("THEIL-SEN ROBUST TREND LINES\n")
cat("============================================\n\n")

cat(
  "|MRBDmax| vs CMCinc:",
  "intercept =",
  coef(ts_MRBD)[1],
  ", slope =",
  coef(ts_MRBD)[2],
  "\n"
)

cat(
  "SRBSmax vs CMCinc:",
  "intercept =",
  coef(ts_SRBS)[1],
  ", slope =",
  coef(ts_SRBS)[2],
  "\n"
)

cat(
  "DeltaSRBS-MRBD vs CMCinc:",
  "intercept =",
  coef(ts_DELTA)[1],
  ", slope =",
  coef(ts_DELTA)[2],
  "\n"
)


## ============================================================
## 9) Scatter plots + Theil-Sen trend lines
## ============================================================


## ------------------------------------------------------------
## 9.1 |MRBDmax| vs CMCinc
## ------------------------------------------------------------

p_MRBD <- ggplot(
  df1,
  aes(
    x = MRBD_mag,
    y = CMCinc
  )
) +
  
  geom_point(
    size = 2
  ) +
  
  geom_abline(
    intercept = coef(ts_MRBD)[1],
    slope = coef(ts_MRBD)[2],
    linewidth = 0.8
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  labs(
    title = "|MRBDmax| vs CMCinc",
    
    subtitle = paste0(
      "Spearman rho = ",
      round(sp_MRBD$rho, 3),
      ", 95% CI [",
      round(sp_MRBD$CI_low, 3),
      ", ",
      round(sp_MRBD$CI_high, 3),
      "], permutation p = ",
      format.pval(
        sp_MRBD$p_perm,
        digits = 3,
        eps = 1 / (R_PERM + 1)
      )
    ),
    
    x = "|MRBDmax|",
    y = "CMCinc"
  )


## ------------------------------------------------------------
## 9.2 SRBSmax vs CMCinc
## ------------------------------------------------------------

p_SRBS <- ggplot(
  df1,
  aes(
    x = SRBSmax,
    y = CMCinc
  )
) +
  
  geom_point(
    size = 2
  ) +
  
  geom_abline(
    intercept = coef(ts_SRBS)[1],
    slope = coef(ts_SRBS)[2],
    linewidth = 0.8
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  labs(
    title = "SRBSmax vs CMCinc",
    
    subtitle = paste0(
      "Spearman rho = ",
      round(sp_SRBS$rho, 3),
      ", 95% CI [",
      round(sp_SRBS$CI_low, 3),
      ", ",
      round(sp_SRBS$CI_high, 3),
      "], permutation p = ",
      format.pval(
        sp_SRBS$p_perm,
        digits = 3,
        eps = 1 / (R_PERM + 1)
      )
    ),
    
    x = "SRBSmax",
    y = "CMCinc"
  )


## ------------------------------------------------------------
## 9.3 DeltaSRBS-MRBD vs CMCinc
## ------------------------------------------------------------

p_DELTA <- ggplot(
  df1,
  aes(
    x = Delta_SRBS_MRBD,
    y = CMCinc
  )
) +
  
  geom_point(
    size = 2
  ) +
  
  geom_abline(
    intercept = coef(ts_DELTA)[1],
    slope = coef(ts_DELTA)[2],
    linewidth = 0.8
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  labs(
    title = expression(
      Delta*"SRBS-MRBD vs CMCinc"
    ),
    
    subtitle = paste0(
      "Spearman rho = ",
      round(sp_DELTA$rho, 3),
      ", 95% CI [",
      round(sp_DELTA$CI_low, 3),
      ", ",
      round(sp_DELTA$CI_high, 3),
      "], permutation p = ",
      format.pval(
        sp_DELTA$p_perm,
        digits = 3,
        eps = 1 / (R_PERM + 1)
      )
    ),
    
    x = expression(
      Delta*"SRBS-MRBD = SRBSmax + |MRBDmax|"
    ),
    
    y = "CMCinc"
  )


print(p_MRBD)
print(p_SRBS)
print(p_DELTA)


## ============================================================
## 10) Key results for reporting
## ============================================================

cat("\n============================================\n")
cat("KEY RESULTS FOR REPORTING\n")
cat("============================================\n")


## ΔSRBS-MRBD vs CMCinc
cat(
  "\nDeltaSRBS-MRBD vs CMCinc:\n",
  "Spearman rho = ",
  round(sp_DELTA$rho, 3),
  "\n95% bootstrap CI = [",
  round(sp_DELTA$CI_low, 3),
  ", ",
  round(sp_DELTA$CI_high, 3),
  "]",
  "\nPermutation p = ",
  format.pval(
    sp_DELTA$p_perm,
    digits = 4,
    eps = 1 / (R_PERM + 1)
  ),
  "\n",
  sep = ""
)


## |MRBDmax| vs CMCinc
cat(
  "\n|MRBDmax| vs CMCinc:\n",
  "Spearman rho = ",
  round(sp_MRBD$rho, 3),
  "\n95% bootstrap CI = [",
  round(sp_MRBD$CI_low, 3),
  ", ",
  round(sp_MRBD$CI_high, 3),
  "]",
  "\nPermutation p = ",
  format.pval(
    sp_MRBD$p_perm,
    digits = 4,
    eps = 1 / (R_PERM + 1)
  ),
  "\n",
  sep = ""
)


## SRBSmax vs CMCinc
cat(
  "\nSRBSmax vs CMCinc:\n",
  "Spearman rho = ",
  round(sp_SRBS$rho, 3),
  "\n95% bootstrap CI = [",
  round(sp_SRBS$CI_low, 3),
  ", ",
  round(sp_SRBS$CI_high, 3),
  "]",
  "\nPermutation p = ",
  format.pval(
    sp_SRBS$p_perm,
    digits = 4,
    eps = 1 / (R_PERM + 1)
  ),
  "\n",
  sep = ""
)


## partial |MRBDmax|
cat(
  "\n|MRBDmax| vs CMCinc controlling SRBSmax:\n",
  "Partial Spearman rho = ",
  round(ps_MRBD$rho, 3),
  "\nPermutation p = ",
  format.pval(
    ps_MRBD$p_perm,
    digits = 4,
    eps = 1 / (R_PERM + 1)
  ),
  "\n",
  sep = ""
)


## partial SRBSmax
cat(
  "\nSRBSmax vs CMCinc controlling |MRBDmax|:\n",
  "Partial Spearman rho = ",
  round(ps_SRBS$rho, 3),
  "\nPermutation p = ",
  format.pval(
    ps_SRBS$p_perm,
    digits = 4,
    eps = 1 / (R_PERM + 1)
  ),
  "\n",
  sep = ""
)