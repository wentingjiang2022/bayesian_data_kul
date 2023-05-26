library(HDInterval)
library(rjags)
library(coda)
library(readr)
library(runjags)


#### 2.3 simple model (Question 1-2) ####
ophtalmology <- read_csv("ophtalmology.txt")
ophtalmology <- ophtalmology[2]    # Drop first column as it is only an index column

# Prepare data:
my.data <- list(lutein_2003 = c(6.4, 7.5, 8.4, 9.6, 12.0, 4.2, 3.1, 6.3, 4.4),
                lutein_1999 = c(3.5, 2.9, 4.1, 5.1, 6.4, 1.9, 1.3, 4.1, 2.3),
                N = 9)

# Initial parameters
my.inits <- list(
  list(beta0 = 0.1, beta1 = 0.025, tau=1/0.05),
  list(beta0 = 0.9, beta1 = 0.6, tau = 1/0.3)
)

# Model specification
cat("model
  {
    for (i in 1:N){
      lutein_2003[i] ~ dnorm(mu[i],tau)
      mu[i] <- beta0 + beta1*lutein_1999[i]
    }
    sigma2 <- 1/tau
    beta0 ~ dnorm(0,1.0E-6)
    beta1 ~ dnorm(0,1.0E-6)
    tau ~ dgamma(1.0E-3,1.0E-3)
  }", file="ophtalmology_model.txt")

# Parameters to monitor
parameters <- c("beta0", "beta1", "sigma2")

# Running JAGS:
jags <- jags.model(file ="ophtalmology_model.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 2)

ophtal.sim <- coda.samples(model = jags,
                           variable.names = parameters,
                           n.iter = 5000, 
                           thin = 1)

# Convert osteo.sim into mcmc.list for processing with CODA package
ophtal.mcmc <- as.mcmc.list(ophtal.sim)

# Produce general summary of obtained MCMC sampling
summary(ophtal.mcmc)
ophtal.combined <- combine.mcmc(ophtal.mcmc)
HPDinterval(ophtal.combined)

# Traceplot obtained from CODA functions
par(mfrow=c(3,1))
traceplot(ophtal.mcmc[,1], main = 'Trace of beta0', cex.main = 1.5)
traceplot(ophtal.mcmc[,2], main = 'Trace of beta1', cex.main = 1.5)
traceplot(ophtal.mcmc[,3], main = 'Trace of sigma2', cex.main = 1.5)
# Too many samples on a too small plot --> "zoom in" somewhere at the end of the chain
traceplot(ophtal.mcmc[,1], main = '"Zoomed in" trace of beta0', cex.main = 1.5, xlim = c(4000, 4500))
traceplot(ophtal.mcmc[,2], main = '"Zoomed in" trace of beta1', cex.main = 1.5, xlim = c(4000, 4500))
traceplot(ophtal.mcmc[,3], main = '"Zoomed in" trace of sigma2', cex.main = 1.5, xlim = c(4000, 4500))

# Autocorrelation plot
par(mfrow=c(1,1))
acfplot(ophtal.mcmc, lag.max = 30, ylim = c(-0.2, 1))
# crosscorr.plot(ophtal.mcmc)

# Gelman-Rubin convergence diagnostics
gelman.diag(ophtal.mcmc)
gelman.plot(ophtal.mcmc[,1], main = 'beta0', ask = FALSE)
gelman.plot(ophtal.mcmc[,2], main = 'beta1', ask = FALSE)
gelman.plot(ophtal.mcmc[,3], main = 'sigma2', ask = FALSE)


#### 2.3 Longer model with thinning, burnin (Question 2) ####
# Only barely mention this in the Assignment --> what's important is the next model
cat("model
  {
    for (i in 1:N){
      lutein_2003[i] ~ dnorm(mu[i],tau)
      mu[i] <- beta0 + beta1*(lutein_1999[i])
    }
    sigma2 <- 1/tau
    beta0 ~ dnorm(0,1.0E-6)
    beta1 ~ dnorm(0,1.0E-6)
    tau ~ dgamma(1.0E-3,1.0E-3)
  }", file="ophtalmology_model_semicomp.txt")


# Parameters to monitor
parameters <- c("beta0", "beta1", "sigma2")

# Running JAGS:
jags <- jags.model(file ="ophtalmology_model_semicomp.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 2)
update(jags,2000)
ophtal.sim <- coda.samples(model = jags,
                           variable.names = parameters,
                           n.iter = 90000, 
                           thin = 10)

# Convert osteo.sim into mcmc.list for processing with CODA package
ophtal.mcmc_semicomp <- as.mcmc.list(ophtal.sim)

# Produce general summary of obtained MCMC sampling
summary(ophtal.mcmc_semicomp)
ophtal.combined <- combine.mcmc(ophtal.mcmc_semicomp)
HPDinterval(ophtal.combined)

# Traceplot obtained from CODA functions
par(mfrow=c(3,1))
traceplot(ophtal.mcmc_semicomp)
densplot(ophtal.mcmc_semicomp)

# Autocorrelation plot
par(mfrow=c(1,1))
acfplot(ophtal.mcmc_semicomp)
# crosscorr.plot(ophtal.mcmc_semicomp)

# Gelman-Rubin convergence diagnostics
gelman.diag(ophtal.mcmc_semicomp)
gelman.plot(ophtal.mcmc_semicomp,ask=FALSE)


#### 2.3 Longer model with thinning, burnin, centering ####
cat("model
  {
    for (i in 1:N){
      lutein_2003[i] ~ dnorm(mu[i],tau)
      mu[i] <- beta0 + beta1*(lutein_1999[i]-mean(lutein_1999[]))
    }
    sigma2 <- 1/tau
    beta0 ~ dnorm(0,1.0E-6)
    beta1 ~ dnorm(0,1.0E-6)
    tau ~ dgamma(1.0E-3,1.0E-3)
  }", file="ophtalmology_model_complex.txt")


# Parameters to monitor
parameters <- c("beta0", "beta1", "sigma2")

# Running JAGS:
jags <- jags.model(file ="ophtalmology_model_complex.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 2)
update(jags,2000)
ophtal.sim <- coda.samples(model = jags,
                           variable.names = parameters,
                           n.iter = 90000, 
                           thin = 10)

# Convert osteo.sim into mcmc.list for processing with CODA package
ophtal.mcmc_complex <- as.mcmc.list(ophtal.sim)

# Produce general summary of obtained MCMC sampling
summary(ophtal.mcmc_complex)
ophtal.combined <- combine.mcmc(ophtal.mcmc_complex)
HPDinterval(ophtal.combined)

# Traceplot obtained from CODA functions
par(mfrow=c(3,1))
traceplot(ophtal.mcmc_complex[,1], main = 'Trace of beta0', cex.main = 1.5)
traceplot(ophtal.mcmc_complex[,2], main = 'Trace of beta1', cex.main = 1.5)
traceplot(ophtal.mcmc_complex[,3], main = 'Trace of sigma2', cex.main = 1.5)
# Too many samples on a too small plot --> "zoom in" somewhere at the end of the chain
traceplot(ophtal.mcmc_complex[,1], main = '"Zoomed in" trace of beta0', cex.main = 1.5, xlim = c(60000, 65000))
traceplot(ophtal.mcmc_complex[,2], main = '"Zoomed in" trace of beta1', cex.main = 1.5, xlim = c(60000, 65000))
traceplot(ophtal.mcmc_complex[,3], main = '"Zoomed in" trace of sigma2', cex.main = 1.5, xlim = c(60000, 65000))

# Autocorrelation plot
par(mfrow=c(1,1))
acfplot(ophtal.mcmc_complex)
# crosscorr.plot(ophtal.mcmc_complex)

# Gelman-Rubin convergence diagnostics
gelman.diag(ophtal.mcmc_complex)
gelman.plot(ophtal.mcmc_complex,ask=FALSE)


#### 2.3 Complex Model with prediction (Question 3) ####
cat("model
  {
    for (i in 1:N){
      lutein_2003[i] ~ dnorm(mu[i],tau)
      mu[i] <- beta0 + beta1*(lutein_1999[i]-mean(lutein_1999[]))
    }
    sigma2 <- 1/tau
    beta0 ~ dnorm(0,1.0E-6)
    beta1 ~ dnorm(0,1.0E-6)
    tau ~ dgamma(1.0E-3,1.0E-3)
    
    lutein_2003_new ~ dnorm(mu_new, tau)
    mu_new <- beta0+beta1*(5-mean(lutein_1999[]))
  }", file="ophtalmology_model_complex_pred.txt")


# Parameters to monitor
parameters <- c("lutein_2003_new")

# Running JAGS:
jags <- jags.model(file ="ophtalmology_model_complex_pred.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 2)
update(jags,2000)
ophtal.sim <- coda.samples(model = jags,
                           variable.names = parameters,
                           n.iter = 90000, 
                           thin = 10)

# Convert osteo.sim into mcmc.list for processing with CODA package
ophtal.mcmc_complex_pred <- as.mcmc.list(ophtal.sim)

# Produce general summary of obtained MCMC sampling
summary(ophtal.mcmc_complex_pred)
ophtal.combined <- combine.mcmc(ophtal.mcmc_complex_pred)
HPDinterval(ophtal.combined)


#### 2.3 Predicting for 100 participants (Question 4) ####
cat("model
  {
    for (i in 1:N){
      lutein_2003[i] ~ dnorm(mu[i],tau)
      mu[i] <- beta0 + beta1*(lutein_1999[i]-mean(lutein_1999[]))
    }
    sigma2 <- 1/tau
    beta0 ~ dnorm(0,1.0E-6)
    beta1 ~ dnorm(0,1.0E-6)
    tau ~ dgamma(1.0E-3,1.0E-3)
    
    # Prediction for new 1999 values
    for (j in 1:N_new) {
      y_new[j] ~ dnorm(mu_new[j], tau)
      mu_new[j] <- beta0 + beta1*(x_new[j]-mean(lutein_1999[]))
    }
  }", file="ophtalmology_model_complex_pred_new_data.txt")

# Define new data
N_new <- 100
x_new <- ophtalmology$x

# Combine new data with original data
my.data <- list(lutein_2003 = c(6.4, 7.5, 8.4, 9.6, 12.0, 4.2, 3.1, 6.3, 4.4),
                lutein_1999 = c(3.5, 2.9, 4.1, 5.1, 6.4, 1.9, 1.3, 4.1, 2.3),
                N = 9,
                y_new = rep(NA, N_new),
                x_new = x_new,
                N_new = N_new)

# Parameters to monitor
parameters <- c("beta0", "beta1", "sigma2", "y_new")

# Running JAGS:
jags <- jags.model(file ="ophtalmology_model_complex_pred_new_data.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 2)
update(jags,2000)
ophtal.sim <- coda.samples(model = jags,
                           variable.names = parameters,
                           n.iter = 90000, 
                           thin = 10)

ophtal.mcmc_complex_pred_new_data <- as.mcmc.list(ophtal.sim)

# Produce general summary of obtained MCMC sampling
summary(ophtal.mcmc_complex_pred_new_data)

# Extract posterior samples
sample_chain1 <- as.data.frame(ophtal.sim[[1]])
sample_chain2 <- as.data.frame(ophtal.sim[[2]])
posterior <- rbind(sample_chain1, sample_chain2)

# Initialize vector to store posterior probabilities
posterior_probs <- numeric(N_new)

# Calculate posterior probability for each new observation
for (j in 1:N_new) {
  posterior_probs[j] <- mean(posterior[, 3 + j] > 5, na.rm = TRUE)
}

# Output posterior probabilities
posterior_probs

# Histogram
hist(posterior_probs)

# b)
# Number of participants who have a high posterior probability that they took serum-lutein capsules
sum(posterior_probs > 0.8)

# Calculate posterior probability for each new observation
for (j in 1:N_new) {
  posterior_probs[j] <- mean(posterior[, 3 + j] <= 5, na.rm = TRUE)
}

# Output posterior probabilities
posterior_probs

# Number of participants who have a high posterior probability that they took no serum-lutein capsules
sum(posterior_probs > 0.8)
