#' @title Escalation and De-escalation threshold Generator
#'
#' @description The function `get.rules` generates escalation and de-escalation rules
#' for several interval designs including the CCD, mTPI,
#' BOIN, Keyboard and UMPBI designs.
#' The function was adapted based on the work of Ruitao Lin
#'
#' @param target Numeric: The target dose-limiting toxicity (DLT) rate; default=0.3.
#' @param ncohort Numeric: The total number of cohorts; scalar.
#' @param cohortsize Integer: The number of subjects per cohort; scalar, default=3.
#' @param design Integer: The interval designs of dose finding to be used; scalar. \cr
#' Default values for parameters in designs are used. \cr
#' \itemize{
#'   \item 1 - Central Composite (CCD) Design
#'   \item 2 - Modified Toxicity Probability Interval (mTPI) Design
#'   \item 3 - Bayesian Optimal Interva (BOIN) Design
#'   \item 4 - Keyboard or mTPI2 deisgn
#'   \item 5 - Uniformly Most Powerful Bayesian Interval (UMPBI) Design
#' }
#' @param cutoff.eli Numeric: Cutoff probability of early termination,
#' avoiding the overly toxic dose for safety monitoring; scalar, default=0.95
#'
#' @returns A table/matrix of escalation, de-escalation, and elimination rules for each cohort.
#'
#' @export
#'
#' @examples get.rules(target=0.3,ncohort=10, cohortsize=3, design=2, cutoff.eli=0.95) # this should match Fig 2 of mTPI original paper
get.rules <- function(target, ncohort, cohortsize=3, design=7, cutoff.eli=0.95){
  # Assert design choices input
  valid_designs <- 1:5
  if (!(design %in% valid_designs)) {
    stop("Invalid design choice! Choose from: 1 (CCD), 2 (mTPI), 3 (BOIN), 4 (Keyboard), 5 (UMPBI).")
  }


  # Parameter setting of design
  if (design == 1) {  # CCD
    if (target < 0.3) { lambda1 <- target - 0.09; lambda2 <- target + 0.09 }
    else if (target < 0.4) { lambda1 <- target - 0.1; lambda2 <- target + 0.1 }
    else if (target < 0.45) { lambda1 <- target - 0.12; lambda2 <- target + 0.12 }
    else { lambda1 <- target - 0.13; lambda2 <- target + 0.13 }
  }


  if (design == 2) {# Bayesian mTPI via Beta-binomial Conjugacy
    # prior distribution of p ~ Beta(alpha, beta) of Bin(n,p)
    prior_alpha <- 1
    prior_beta <- 1
    esp <- 0.05
    p.saf <- target - esp
    p.tox <- target + esp

    # Probability Mass function for left tail, middle, and right tail.
    get.pm1 <- function(n, m){pbeta(p.saf,m+1,n-m+1)}
    get.pm2 <- function(n, m){pbeta(p.tox, m+1, n-m+1)-pbeta(p.saf,m+1,n-m+1)}
    get.pm3 <- function(n, m){1- pbeta(p.tox, m+1, n-m+1)}
    # exclusion prob.
    get.ep <- function(n, m){1-pbeta(target, m+1, n-m+1)}
  }


  if (design==3) { # BOIN
    p.saf <- target * 0.6
    p.tox <- target * 1.4
    lambda1 <- log((1-p.saf)/(1-target)) / log(target*(1-p.saf)/(p.saf*(1-target)))
    lambda2 <- log((1-target)/(1-p.tox)) / log(p.tox*(1-target)/(target*(1-p.tox)))
  }

  if (design ==4) {# Keyboard
    esp <- 0.05
    # Calculate probability mass for the acceptable range [m1,m2]
    get.pm1 <- function(p, n, m1){1 - pbinom(m1, n, p)}
    get.pm2 <- function(p, n, m1, m2){pbinom(m1, n, p) + 1 - pbinom(m2 - 1, n, p)}
    get.pm3 <- function(p, n, m2){pbinom(m2 - 1, n, p)}

    p.saf <- target - esp
    p.tox <- target + esp
  }
  if (design==5) {c<-log(1.1)/3} # UMPBI constant

  # Numerical search for decision boundaries
  npts = ncohort*cohortsize;
  ntrt = NULL; b.e = NULL; b.d = NULL; elim = NULL;

  for(n in (1:ncohort)*cohortsize) {
    # total enrollment over time
    ntrt = c(ntrt, n)
    if (design %in% c(1, 3)) { # CCD or BOIN
      cutoff1 = floor(n*lambda1)
      cutoff2 = floor(n*lambda2)+1
    }

    if (design == 2){ # mTPI
      cutoff1=-1; cutoff2=n+1
      for (i in 0:n){
        upm1 <- get.pm1(n,i) / p.saf
        upm2 <- get.pm2(n,i) / (p.tox - p.saf)
        upm3 <- get.pm3(n,i) / (1 - p.tox)
        if (upm1 == max(c(upm1,upm2,upm3))){cutoff1 <- max(i,cutoff1)}
        if (upm3 == max(c(upm1,upm2,upm3))){cutoff2 <- min(i,cutoff2)}
        }
      }


    if (design == 4) { # Keyboard
      error.min=3;
      for (m1 in 0:floor(target*n)){
        for (m2 in ceiling(target*n):n){
          if (design==2) {
            print(p.saf)
            error1 = integrate(get.pm1, lower=p.saf, upper=p.tox, n, m1, m2)$value/(p.tox-p.saf);
            error2 = integrate(get.pm2, lower=0, upper=p.saf, n, m1)$value/p.saf;
            error3 = integrate(get.pm3, lower=p.tox, upper=1, n, m2)$value/(1-p.tox);
          }
          if (design==4) {
            epsilon=p.tox-p.saf
            error1 = integrate(get.pm1, lower=p.saf, upper=p.tox, n, m1, m2)$value/epsilon;
            error2 = integrate(get.pm2, lower=p.saf-epsilon, upper=p.saf, n, m1)$value/epsilon;
            error3 = integrate(get.pm3, lower=p.tox, upper=p.tox+epsilon, n, m2)$value/epsilon;
          }
          error = error1 + error2 + error3
          if(error<error.min) {
            error.min = error
            cutoff1 = m1
            cutoff2 = m2
          }
        }
      }
    }

    if (design==5) { # UMPBI
      gamma = exp(c*sqrt(n))
      df <- function(p,gamma,n,target){
        l <- (n/(1-p)) / (log(p/(1-p))-log(target/(1-target)))
        l <- l- (log(gamma)-n*log((1-p)/(1-target)))/p/(1-p) / (log(p/(1-p))-log(target/(1-target)))^2
      }
      p.tox = uniroot(df,c(target+0.0001,target*2),gamma=gamma,n=n,target=target)$root
      p.saf = uniroot(df,c(0.0001,target-0.0001),gamma=gamma,n=n,target=target)$root
      lambda1  = (log((1-p.saf)/(1-target))-log(gamma)/n) / log(target*(1-p.saf)/(p.saf*(1-target)));
      lambda2  = (log((1-target)/(1-p.tox))+log(gamma)/n) / log(p.tox*(1-target)/(target*(1-p.tox)));
      cutoff1 = floor(n*lambda1)
      cutoff2 = floor(n*lambda2)+1
    }

    # Store escalation and deescalation rules
    b.e = c(b.e, cutoff1);
    b.d = c(b.d, cutoff2);

    # Compute elimination boundaries
    if (n < 3) {
    } else {
      elim_threshold <- NA
      for (ntox in 3:n) {
        if (1 - pbeta(target, ntox + 1, n - ntox + 1) > cutoff.eli) {
          elim_threshold <- ntox
          break
        }
      }
      elim <- c(elim, elim_threshold)
    }
  }
  # Adjust de-escalation rules to not exceed elimination boundaries
  b.d <- ifelse(!is.na(elim) & b.d > elim, elim, b.d)


  # Construct output matrix
  rules <- rbind(ntrt, b.e, b.d, elim)
  rownames(rules) <- c("Number of patients treated",
                       "Escalate   if # of DLT <=",
                        "Deescalate if # of DLT >=",
                        "Terminate  if # of DLT >=")
  colnames(rules) <- rep("", ncohort)

  return(rules);
}

#> get.rules(target=0.3,ncohort=10, cohortsize=3, design=2, cutoff.eli=0.95)
# Number of patients treated 3 6 9 12 15 18 21 24 27 30
# Escalate   if # of DLT <=  0 1 1  2  2  3  4  4  5  6
# Deescalate if # of DLT >=  2 4 5  6  8  9 10 11 12 14
# Terminate  if # of DLT >=  3 4 5  7  8  9 10 11 12 14
