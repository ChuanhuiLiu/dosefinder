#' @title Pool Adjacent Violator Algorithm (PAVA)
#' @description an efficient isotonic regression/smoothing method that enforces a monotonic constraint on a sequence of values, e.g., toxicity probabilities for dose-response modeling. The function was adapted based on the work of Ruitao Lin
#' @param x Numeric vector: Input values.
#' @param wt Numeric vector: Weights (default is uniform weight).
#' @return Numeric vector: Isotonically adjusted values.
#' @examples # toxicity_rates <- c(0.05, 0.10, 0.15, 0.14, 0.20) \cr
#' # pava(toxicity_rates) \cr
#' #[1] 0.05 0.10 0.15 0.15 0.20 \cr
pava <- function(x, wt = rep(1, length(x))) {
  n <- length(x)
  if (n <= 1) return(x)
  if (any(is.na(x)) || any(is.na(wt))) stop("Missing values in 'x' or 'wt' not allowed")

  lvlsets <- seq_len(n)
  repeat {
    # Identify violations of monotonicity
    viol <- diff(x) < 0
    if (!any(viol)) break
    # Find the first violation
    i <- min(which(viol))

    # Pool adjacent violators
    ilvl <- lvlsets %in% c(lvlsets[i], lvlsets[i + 1])
    x[ilvl] <- sum(x[ilvl] * wt[ilvl]) / sum(wt[ilvl])  # Adjust values
    lvlsets[ilvl] <- lvlsets[i]
  }
  return(x)
}
