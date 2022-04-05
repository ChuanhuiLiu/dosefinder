#' @title Operating Characteristics for Bayesian Dose-Finding Designs
#' @description Simulates clinical trials using interval-based dose-finding method. The function was adapted based on the work of Ruitao Lin
#' @param target Numeric: Target toxicity rate.
#' @param p.true Numeric vector: True toxicity probability for each dose.
#' @param ncohort Integer: Number of patient cohorts.
#' @param cohortsize Integer: Number of patients per cohort.
#' @param startdose Integer: Initial dose level.
#' @param design Integer: Design choice (1=CCD, 2=mTPI, 3=BOIN, 4=Keyboard, 5=UMPBI).
#' @param design Integer: Estimation choice (1=PAVA, 2=BLRM).
#' @param cutoff.eli Numeric: Cutoff probability for dose elimination (default = 0.95).
#' @param ntrial Integer: Number of simulated trials.
#' @return A data frame summarizing trial results; A list of summary
#' @export
#' @examples p.true<-c(0.08,0.10,0.20,0.30,0.45,0.60); \br
#' summary.simulation(target=0.3,p.true=p.true,ncohort=12,cohortsize=3,type =2,design=5,cutoff.eli=0.95)
summary.simulation<-function(target, p.true, ncohort, cohortsize, startdose = 1, design = 5, type =2, cutoff.eli = 0.95, ntrial = 1000) {
  # Assert design choices input
  valid_designs <- 1:7
  if (!(design %in% valid_designs)) {
    stop("Invalid design choice! Choose from: 1 (CCD), 2 (mTPI), 3 (BOIN), 4 (Keyboard), 5 (UMPBI), 6 (CRM), 7 (Hybrid).")
  }

  # Set random seed for reproducibility
  set.seed(6)

  # Number of dose levels
  ndose <- length(p.true)
  npts <- ncohort * cohortsize

  # Matrices to store simulation results
  Y <- matrix(0, nrow = ntrial, ncol = ndose)  # Toxicity counts
  N <- matrix(0, nrow = ntrial, ncol = ndose)  # Patient counts
  dselect <- numeric(ntrial)  # Selected MTD per trial

  # Obtain dose escalation and de-escalation rules
  boundaries <- get.rules(target, ncohort, cohortsize, design, cutoff.eli)
  b.e <- boundaries[2,]   # Escalation boundary
  b.d <- boundaries[3,]   # De-escalation boundary
  b.elim <- boundaries[4,]  # Elimination boundary

  # Simulate trials
  for (trial in seq_len(ntrial)) {
    y <- numeric(ndose)  # Toxicity counts per dose
    n <- numeric(ndose)  # Patients treated per dose
    d <- startdose  # Starting dose level
    earlystop <- FALSE  # Track early termination
    elimi <- rep(0, ndose)  # Dose elimination flag

    for (i in seq_len(ncohort)) {

      # Simulate toxicity outcome
      y[d] <- y[d] + sum(runif(cohortsize) < p.true[d])
      n[d] <- n[d] + cohortsize
      nc <- n[d] / cohortsize  # Normalize patient count

      # Check for dose elimination
      if (!is.na(b.elim[nc]) && y[d] >= b.elim[nc]) {
        elimi[d:ndose] <- 1
        if (d == 1) { earlystop <- TRUE; break }
      }

      # Dose escalation/de-escalation decision
      if (y[d] <= b.e[nc] && d != ndose && elimi[d + 1] == 0) {
        d <- d + 1  # Escalate
      } else if (y[d] >= b.d[nc] && d != 1) {
        d <- d - 1  # De-escalate
      }
    }

    # Store trial results
    Y[trial, ] <- y
    N[trial, ] <- n
    a <-if (earlystop) 99 else select.mtd(target, y, n, type, cutoff.eli)
    dselect[trial] <- a
  }

  # Compute summary statistics
  selpercent <- colMeans(dselect == matrix(1:ndose, nrow = ntrial, ncol = ndose, byrow = TRUE)) * 100
  patient_counts <- colMeans(N)
  toxicity_counts <- colMeans(Y)

  summary_text <- c(
    paste("Total average toxicities:", formatC(sum(Y) / ntrial, digits = 1, format = "f")),
    paste("Total average patients:", formatC(sum(N) / ntrial, digits = 1, format = "f"))
  )

  return(list(Selec_Percent = formatC(selpercent, digits = 1, format = "f"),
              Avg_Patients = formatC(patient_counts, digits = 1, format = "f"),
              Avg_Toxicities = formatC(toxicity_counts, digits = 1, format = "f"),
              Summary = summary_text))
}

#p.true<-c(0.08,0.10,0.20,0.30,0.45,0.60)
#summary.simulation(target=0.3,p.true=p.true,ncohort=12,cohortsize=3,design=2,type=1, cutoff.eli=0.95, ntrial = 1000)
# ~37% prob. of selecting the lower dose; ONLY 46% prob. of finding the correct one
 summary.simulation(target=0.3,p.true=p.true,ncohort=12,cohortsize=3,design=2,type=2, cutoff.eli=0.95, ntrial = 20) # less ntrial due to MCMC involved
# ~35% prob. of selecting the lower dose; 60% prob. of finding the correct one (14% improvement)
