#' @title Bayesian Logistic Regression Model (BLRM) for modeling the posterior toxicity probability
#' @description The function `estimate.blrm` estimates Posterior Toxicity Probabilities using observed dose-limiting toxicities (DLT) y and number of enrolled patients n.
#' The function was adapted based on the work of Ruitao Lin.
#' @param target Numeric: Target toxicity rate.
#' @param y Integer vector: Number of accumulated dose-limiting toxicities (DLTs) at each dose.
#' @param n Integer vector: Number of accumulated patients treated at each dose.
#' @param d Numeric vector: vector of dose level
#' @param plot Boolean: show plot of the estimated posterior probability
#' @return a dataframe containing Posterior Mean and 95% Credible Interval for each dose
#' @export
#' @examples dose <- c(1, 2, 3, 4, 5);y <- c(1, 2, 3, 4, 5);n <- c(3, 6, 9, 12, 15); estimate.blrm(y,n,dose)
#' library(rstanarm)  # uncomment this for standalone use
#' library(ggplot2)
estimate.blrm <- function(y,n,d,plot=FALSE){

  # Normalize doses for stability (log transformation)
  dose_scaled <- log(d / min(d))
  # Bayesian Logistic Regression Model (BLRM)
  suppressWarnings(capture.output(
  blrm_model <- stan_glm(
    formula = cbind(y, n - y) ~ dose_scaled,  # Binomial model for toxicity probability
    family = binomial(link = "logit"),  # Logistic regression
    prior = normal(0, 2.5, autoscale = TRUE),  # Weakly informative prior
    prior_intercept = normal(0, 5, autoscale = TRUE),  # Prior on intercept
    iter = 4000, chains = length(n), seed = 1234,refresh = -1 # MCMC settings
  )))

  # Compute Posterior Mean and 95% Credible Interval for Each Dose
  posterior_toxicity_prob <- posterior_epred(blrm_model)  # Matrix of posterior probabilities

  # Compute posterior mean per dose
  posterior_means <- colMeans(posterior_toxicity_prob)

  # Compute 95% credible intervals properly (matching doses)
  cred_intervals <- apply(posterior_toxicity_prob, 2, quantile, probs = c(0.025, 0.975))


  # Plot Posterior Toxicity Probability vs. Dose
  if (plot == TRUE){
    results <- data.frame(
      Dose = d,
      Observed_DLT = y,
      Enrolled = n,
      Posterior_Mean_Toxicity = colMeans(posterior_toxicity_prob),
      Lower_95_CI = cred_intervals[,1],
      Upper_95_CI = cred_intervals[,2]
    )
    ggplot(results, aes(x = Dose, y = Posterior_Mean_Toxicity)) +
      geom_point(size = 3, color = "blue") +
      geom_errorbar(aes(ymin = Lower_95_CI, ymax = Upper_95_CI), width = 0.2) +
      labs(title = "Bayesian Posterior Toxicity Probability (BLRM)",
           x = "Dose Level",
           y = "Estimated Toxicity Probability") +
      theme_minimal()
  }
  return(list(Dose = d,Observed_DLT = y,Enrolled = n,
              Posterior_Mean_Toxicity = posterior_means,
              Lower_95_CI = cred_intervals[1, ],
              Upper_95_CI = cred_intervals[2, ]))
}
