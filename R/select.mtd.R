#' @title Select Maximum Tolerated Dose (MTD) in a Dose-Finding Trial
#' @description The function `select.mtd` estimates Posterior Toxicity Probabilities and determines the MTD based on toxicity data.
#' The function was adapted based on the work of Ruitao Lin.
#' @param target Numeric: Target toxicity rate.
#' @param y Integer vector: Number of accumulated dose-limiting toxicities (DLTs) at each dose.
#' @param n Integer vector: Number of accumulated patients treated at each dose.
#' @param type Integer: The method for estimating posterior toxicity probabilities. \cr
#' \itemize{
#'   \item 1 - PAVA regression with a weakly informative conjugate prior with Beta(0.05,0.05)
#'   \item 2 - Bayesian Logistic Regression Model (BLRM) with a weakly informative normal prior
#' }
#' @param cutoff.eli Numeric: Threshold probability for eliminating an overly toxic dose (default = 0.95).
#' @return Integer: Selected MTD dose level (or 99 if no dose is safe).
#' @export
select.mtd <- function(target, y, n, type=2,cutoff.eli=0.95)
{
  # Step 1: Dose Elimination Based on Excessive Toxicity
  ndose <- length(n)
  elimi <- rep(0, ndose)  # Track eliminated doses

  for (i in seq_len(ndose)) {
    if (n[i] > 2 && (1 - pbeta(target, y[i] + 1, n[i] - y[i] + 1)) > cutoff.eli) {
      elimi[i:ndose] <- 1  # Eliminate dose and all higher doses
      break
    }
  }
  # If the lowest dose is eliminated, no dose should be selected
  if (elimi[1] == 1) return(99)
  # Step 2: Estimate Posterior Toxicity Probabilities
  # Identify the highest admissible (non-eliminated) dose
  nadmis <- min(max(which(elimi == 0)), max(which(n != 0)))
  if (type==1){ # PAVA
  # Compute posterior means and variances using Beta(0.005, 0.005) prior
  phat <- (y[1:nadmis] + 0.005) / (n[1:nadmis] + 0.01)
  phat_var <- (y[1:nadmis] + 0.005) * (n[1:nadmis] - y[1:nadmis] + 0.005) /
    ((n[1:nadmis] + 0.01)^2 * (n[1:nadmis] + 0.02))
  # Apply Isotonic Regression (PAVA)
  phat <- pava(phat, wt = 1 / phat_var)  # Ensure monotonicity
  phat <- phat + seq_len(nadmis) * 1E-10  # Break ties with small increments

  } else if (type == 2){ #BLRM
    phat <- estimate.blrm(y[1:nadmis],n[1:nadmis],1:nadmis)$Posterior_Mean_Toxicity
  }
  # Step 4: Select Dose Closest to Target
  return(which.min(abs(phat - target)))  # Select MTD closest to the target toxicity rate
}

