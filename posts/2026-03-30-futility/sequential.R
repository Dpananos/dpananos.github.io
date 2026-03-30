control_mean <- 1
treatment_mean <- 1

prior_mu <- 0
prior_sd <- 0.05

N <- 100000
confidence_level <- 0.95
sample_size_lower_bound <- 1e4

# =============================================================================
# Ported from ffe: inference_sequential.go
# =============================================================================

# Core bound calculation (normal mixture sequential probability ratio test)
sequential_bound_log_space <- function(observations, negative_2_log_alpha, sample_size_lower_bound) {
  rho <- sample_size_lower_bound / (negative_2_log_alpha + log1p(negative_2_log_alpha)) - 1
  tau <- observations + rho
  sqrt((tau / observations) * (log(tau / rho) + negative_2_log_alpha))
}

# Sequential CI half-width: bound * SE
sequential_ci_bound <- function(n, se, confidence_level, sample_size_lower_bound) {
  neg2la <- -2 * log(1 - confidence_level)
  sequential_bound_log_space(n, neg2la, sample_size_lower_bound) * se
}

# =============================================================================
# Lift + SE helper
# =============================================================================

compute_lift <- function(treatment_yhat, control_yhat, n_per_arm) {
  cofv_c <- 1 / sqrt(n_per_arm) / control_yhat
  cofv_t <- 1 / sqrt(n_per_arm) / treatment_yhat
  lift <- treatment_yhat / control_yhat - 1
  se   <- sqrt((1 + lift)^2 * (cofv_c^2 + cofv_t^2))
  list(lift = lift, se = se)
}

# Is the lift significant under the sequential CI?
is_significant_sequential <- function(lift, se, n, confidence_level, sample_size_lower_bound) {
  bound <- sequential_ci_bound(n, se, confidence_level, sample_size_lower_bound)
  abs(lift) > bound
}

# =============================================================================
# Simulate initial experiment
# =============================================================================

n_per_arm <- N %/% 2

control_yhat <- rnorm(1, control_mean, 1 / sqrt(n_per_arm))
treatment_yhat <- rnorm(1, treatment_mean, 1 / sqrt(n_per_arm))

obs <- compute_lift(treatment_yhat, control_yhat, n_per_arm)
(is_significant <- is_significant_sequential(obs$lift, obs$se, N, confidence_level, sample_size_lower_bound))

# =============================================================================
# Posterior from interim data
# =============================================================================

posterior_V <- 1 / (1 / obs$se^2 + 1 / prior_sd^2)
posterior_mean <- (obs$lift / obs$se^2 + prior_mu / prior_sd^2) * posterior_V

# =============================================================================
# Simulate continuing the experiment
# =============================================================================

Kmore <- 100000
kmore_per_arm <- Kmore %/% 2
combined_per_arm <- n_per_arm + kmore_per_arm
N_combined <- N + Kmore

replicate(1000, {

  new_lift <- rnorm(1, posterior_mean, sqrt(posterior_V))
  control_yhat2 <- rnorm(1, control_mean, 1 / sqrt(kmore_per_arm))
  treatment_yhat2 <- rnorm(1, control_yhat2 * (1 + new_lift), 1 / sqrt(kmore_per_arm))

  new_treatment_yhat <- weighted.mean(c(treatment_yhat, treatment_yhat2), c(n_per_arm, kmore_per_arm))
  new_control_yhat <- weighted.mean(c(control_yhat, control_yhat2), c(n_per_arm, kmore_per_arm))

  combined <- compute_lift(new_treatment_yhat, new_control_yhat, combined_per_arm)
  is_significant_sequential(combined$lift, combined$se, N_combined, confidence_level, sample_size_lower_bound)

}) -> r

mean(r)

# =============================================================================
# Calibration check
# =============================================================================
# For many simulated experiments:
#   1. Draw a true lift from the prior
#   2. Simulate N observations (interim), compute posterior
#   3. Use replicate loop to get predicted P(significant if we continue)
#   4. Simulate the actual Kmore additional data and check if we get significance
#   5. Compare predicted vs observed with rms::val.prob

set.seed(42)
n_outer <- 10000
n_inner <- 1000

results <- t(vapply(seq_len(n_outer), function(i) {

  # 1. Draw true lift from prior
  true_lift <- rnorm(1, prior_mu, prior_sd)

  # 2. Simulate interim data
  c_yhat <- rnorm(1, control_mean, 1 / sqrt(n_per_arm))
  t_yhat <- rnorm(1, control_mean * (1 + true_lift), 1 / sqrt(n_per_arm))

  interim <- compute_lift(t_yhat, c_yhat, n_per_arm)

  # Posterior
  post_V <- 1 / (1 / interim$se^2 + 1 / prior_sd^2)
  post_M <- (interim$lift / interim$se^2 + prior_mu / prior_sd^2) * post_V

  # 3. Predicted P(significant) via inner simulation
  predicted_p <- mean(replicate(n_inner, {
    sim_lift <- rnorm(1, post_M, sqrt(post_V))
    c_yhat2 <- rnorm(1, control_mean, 1 / sqrt(kmore_per_arm))
    t_yhat2 <- rnorm(1, c_yhat2 * (1 + sim_lift), 1 / sqrt(kmore_per_arm))

    sim_t <- weighted.mean(c(t_yhat, t_yhat2), c(n_per_arm, kmore_per_arm))
    sim_c <- weighted.mean(c(c_yhat, c_yhat2), c(n_per_arm, kmore_per_arm))

    sim_combined <- compute_lift(sim_t, sim_c, combined_per_arm)
    is_significant_sequential(sim_combined$lift, sim_combined$se, N_combined, confidence_level, sample_size_lower_bound)
  }))

  # 4. Simulate actual continuation and check significance
  c_yhat_future <- rnorm(1, control_mean, 1 / sqrt(kmore_per_arm))
  t_yhat_future <- rnorm(1, control_mean * (1 + true_lift), 1 / sqrt(kmore_per_arm))

  actual_t <- weighted.mean(c(t_yhat, t_yhat_future), c(n_per_arm, kmore_per_arm))
  actual_c <- weighted.mean(c(c_yhat, c_yhat_future), c(n_per_arm, kmore_per_arm))

  actual_combined <- compute_lift(actual_t, actual_c, combined_per_arm)
  actual_sig <- is_significant_sequential(actual_combined$lift, actual_combined$se, N_combined, confidence_level, sample_size_lower_bound)

  c(predicted = predicted_p, actual = as.numeric(actual_sig))
}, numeric(2)))

rms::val.prob(results[, "predicted"], results[, "actual"])

