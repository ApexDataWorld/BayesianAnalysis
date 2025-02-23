# ====================================================
# ST440/540 Applied Bayesian Analysis – Exam 1
# Due Date: February 17, 2025
# Author: Saurabh Gupta
# ====================================================

# ----- 1. Likelihood and Graphical Assessment -----
library(dplyr)
data("storms")
# Only analyze hurricanes
storms <- storms[!is.na(storms$category),]  

# Extract basic variables from storms
year <- storms$year
name <- paste0(storms$name, "_", storms$year)
wind <- storms$wind

# Create unique storm identifiers and extract max wind speed per storm
uni <- unique(name)
n <- length(uni)
year_storm <- rep(0, n)
mxspd <- rep(0, n)
for(i in 1:length(uni)) {
  # For each unique storm, use all observations with that identifier
  year_storm[i] <- min(year[name == uni[i]])
  mxspd[i] <- max(wind[name == uni[i]])
}

# Compute annual hurricane counts
year_counts <- as.numeric(names(table(year_storm)))
cnt <- as.vector(table(year_storm))

# Plot maximum wind speed by storm year and annual counts (base R)
plot(year_storm, mxspd, xlab="Year", ylab="Max wind speed (knots)",
     main="Max Wind Speed per Storm", pch=16, col="blue")
plot(year_counts, cnt, xlab="Year", ylab="Number of hurricanes",
     main="Annual Hurricane Counts", type="o", col="green", pch=16)


# Split data into Period1: 1975 - 1999 and Period2: 2000 - 2021
period <- ifelse(year_storm < 2000, "Period1", "Period2")
# For annual counts, find period based on year
period_year <- ifelse(year_counts < 2000, "Period1", "Period2")

# For maximum wind speed analysis, find log transformed mxspd
log_mxspd <- log(mxspd)

# Plot histograms and QQ plots for each period (raw and log transformed)
par(mfrow=c(2,2))
for(p in c("Period1", "Period2")) {
  idx <- (period == p)
  raw_data <- mxspd[idx]
  log_data <- log_mxspd[idx]
  
  hist(raw_data, main=paste("Histogram of mxspd", p),
       xlab="Max Wind Speed (knots)", col="blue", breaks=20)
  qqnorm(raw_data, main=paste("QQ Plot of mxspd", p))
  qqline(raw_data, col="red", lwd=2)
  
  hist(log_data, main=paste("Histogram of log(mxspd)", p),
       xlab="log(mxspd)", col="green", breaks=20)
  qqnorm(log_data, main=paste("QQ Plot of log(mxspd)", p))
  qqline(log_data, col="red", lwd=2)
}
par(mfrow=c(1,1))

# ------------------------------------------------------------
# 2. Conjugate Uninformative Prior Distributions
# ------------------------------------------------------------
#
# For Maximum Wind Speed (log-transformed):
#   Likelihood: log(mxspd_i) ~ N(μ, σ²)
#
#   Priors (vague, conjugate):
#     μ | σ² ~ N(μ₀, σ²/τ₀)   where μ₀ = mean(all_log_mxspd), τ₀ = 0.001
#     σ²    ~ Inverse-Gamma(α₀, β₀)   with α₀ = 0.001, β₀ = 0.001
#
# For Hurricane Counts:
#   Likelihood: Y ~ Poisson(λ)
#
#   Prior:
#     λ ~ Gamma(a₀, b₀) with a₀ = 0.01 and b₀ = 0.01

# Set priors for log (mxspd) analysis
all_log <- log_mxspd
mu0 <- mean(all_log)
tau0 <- 0.001
alpha0 <- 0.001
beta0 <- 0.001

# Set priors for hurricane counts
a0_pois <- 0.01
b0_pois <- 0.01

# ------------------------------------------------------------
# 3. Posterior Summaries and Testing for Difference Between Periods
# ------------------------------------------------------------
# For log(mxspd):
#
#   Posterior (conjugate for Normal likelihood):
#     τₙ = τ₀ + n
#     μₙ = (τ₀μ₀ + n·x̄)/(τ₀+n)
#     αₙ = α₀ + n/2
#     βₙ = β₀ + 0.5*(n-1)*s² + (τ₀·n*(x̄-μ₀)²)/(2*(τ₀+n))
#
# Find posterior samples for Period1 and Period2 and then find the difference.
#
# For Hurricane Counts
#   Posterior:
#     λ | Y ~ Gamma(a₀ + ΣY, b₀ + n)
#
# Find posterior samples for each period and find the difference.

# Find posterior parameters for log(mxspd)
calculate_posterior <- function(data, mu0, tau0, alpha0, beta0) {
  n <- length(data)
  xbar <- mean(data)
  s2 <- var(data)
  tau_n <- tau0 + n
  mu_n <- (tau0 * mu0 + n * xbar) / tau_n
  alpha_n <- alpha0 + n / 2
  beta_n <- beta0 + 0.5 * (n-1) * s2 + (tau0 * n * (xbar - mu0)^2) / (2 * tau_n)
  return(list(mu_n=mu_n, tau_n=tau_n, alpha_n=alpha_n, beta_n=beta_n, n=n))
}

# Separate log(mxspd) for each period
log_period1 <- log_mxspd[period == "Period1"]
log_period2 <- log_mxspd[period == "Period2"]

# Compute posterior parameters for each period
posterior_log_period1 <- calculate_posterior(log_period1, mu0,tau0,alpha0,beta0)
posterior_log_period2 <- calculate_posterior(log_period2, mu0,tau0,alpha0,beta0)

set.seed(123)
n_samples <- 10000

# For Period 1
sigma2_period1 <- 1 / rgamma(n_samples, shape = posterior_log_period1$alpha_n, 
                        rate = posterior_log_period1$beta_n)
mu_samples_period1 <- rnorm(n_samples, mean = posterior_log_period1$mu_n, 
                       sd = sqrt(sigma2_period1 / posterior_log_period1$tau_n))

# For Period 2
sigma2_period2 <- 1 / rgamma(n_samples, shape = posterior_log_period2$alpha_n, 
                        rate = posterior_log_period2$beta_n)
mu_samples_period2 <- rnorm(n_samples, mean = posterior_log_period2$mu_n, 
                       sd = sqrt(sigma2_period2 / posterior_log_period2$tau_n))

# Posterior difference (Period2 - Period1)
diff_mu <- mu_samples_period2 - mu_samples_period1
point_estimate_diff <- mean(diff_mu)
standard_error_diff <- sd(diff_mu)
print(paste("Posterior difference in mean log(mxspd):", point_estimate_diff))
print(paste("Standard error:", standard_error_diff))

# For hurricane counts
# Create indices for each period
index1 <- year_counts < 2000
index2 <- year_counts >= 2000

# Separate annual counts by period
count_period1 <- cnt[index1]
count_period2 <- cnt[index2]

# Number of years in each period
n_period1 <- sum(index1)
n_period2 <- sum(index2)

# Compute posterior parameters for λ (Poisson likelihood)
posterior_count_period1 <- list(a_post = a0_pois + sum(count_period1),
                    b_post = b0_pois + n_period1)
posterior_count_period2 <- list(a_post = a0_pois + sum(count_period2),
                    b_post = b0_pois + n_period2)

# Draw posterior samples for λ in each period
lambda_samples_p1 <- rgamma(n_samples, shape = posterior_count_period1$a_post, 
                            rate = posterior_count_period1$b_post)
lambda_samples_p2 <- rgamma(n_samples, shape = posterior_count_period2$a_post, 
                            rate = posterior_count_period2$b_post)

# Posterior difference for hurricane counts:
diff_lambda <- lambda_samples_p2 - lambda_samples_p1
point_estimate_lambda <- mean(diff_lambda)
standard_error_lambda <- sd(diff_lambda)
print(paste("Posterior difference in λ :", point_estimate_lambda))
print(paste("Standard error:", standard_error_lambda))

par(mfrow=c(1,2))
hist(diff_mu, breaks=50, main="Posterior Difference in log(mxspd) 
     (Period2 - Period1)",
     xlab="Difference (log scale)", col="gray")
abline(v=0, col="red", lwd=2)

hist(diff_lambda, breaks=50, main="Posterior Difference in Hurricane Rate (λ)",
     xlab="Difference in λ", col="lightgray")
abline(v=0, col="red", lwd=2)
par(mfrow=c(1,1))

# ------------------------------------------------------------
# 4. Sensitivity Analysis: Check if results are sensitive to priors
# ------------------------------------------------------------
# For log(mxspd): Use alternative hyperparameters
tau0_alternate <- 0.01; 
alpha0_alternate <- 0.01; 
beta0_alternate <- 0.01

posterior_log_period1_alt <- calculate_posterior(log_period1, mu0, 
                             tau0_alternate, alpha0_alternate, beta0_alternate)
posterior_log_period2_alt <- calculate_posterior(log_period2, mu0, 
                             tau0_alternate, alpha0_alternate, beta0_alternate)

# Find posterior samples for μ in each period 
sigma2_period1_alt <- 1 / rgamma(n_samples, 
                            shape = posterior_log_period1_alt$alpha_n, 
                            rate = posterior_log_period1_alt$beta_n)
mu_samples_period1_alt <- rnorm(n_samples, mean =posterior_log_period1_alt$mu_n, 
                sd = sqrt(sigma2_period1_alt / posterior_log_period1_alt$tau_n))

sigma2_period2_alt <- 1 / rgamma(n_samples, 
                            shape = posterior_log_period2_alt$alpha_n, 
                            rate = posterior_log_period2_alt$beta_n)
mu_samples_period2_alt <- rnorm(n_samples, mean =posterior_log_period2_alt$mu_n, 
                sd = sqrt(sigma2_period2_alt / posterior_log_period2_alt$tau_n))

diff_mu_alt <- mu_samples_period2_alt - mu_samples_period1_alt
print(paste("Sensitivity Analysis for log(mxspd) - Difference Alternate Prior:",
            mean(diff_mu_alt)))
print(paste("Sensitivity Analysis for log(mxspd)- Std Error:", sd(diff_mu_alt)))

# For hurricane counts, alternative prior hyperparameters
a0_pois_alt <- 0.1
b0_pois_alt <- 0.1

index1 <- year_counts < 2000
index2 <- year_counts >= 2000
n_period1 <- sum(index1)
n_period2 <- sum(index2)
count_period1 <- cnt[index1]
count_period2 <- cnt[index2]
total_period1 <- sum(count_period1)
total_period2 <- sum(count_period2)

# Find posterior parameters for each period
posterior_count_period1_alt <- list(a_post = a0_pois_alt + total_period1,
                        b_post = b0_pois_alt + n_period1)
posterior_count_period2_alt <- list(a_post = a0_pois_alt + total_period2,
                        b_post = b0_pois_alt + n_period2)

# Posterior samples for lambda
lambda_samples_p1_alt <- rgamma(n_samples, 
                                shape = posterior_count_period1_alt$a_post, 
                                rate = posterior_count_period1_alt$b_post)
lambda_samples_p2_alt <- rgamma(n_samples, 
                                shape = posterior_count_period2_alt$a_post, 
                                rate = posterior_count_period2_alt$b_post)

diff_lambda_alt <- lambda_samples_p2_alt - lambda_samples_p1_alt
print(paste("Sensitivity Analysis (Hurricane Counts) - Difference (Alt Prior):", 
            mean(diff_lambda_alt)))
print(paste("Sensitivity Analysis (Hurricane Counts) - Std Error:", 
            sd(diff_lambda_alt)))

