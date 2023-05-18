library(HDInterval)
library(rjags)
library(coda)
library(readr)
library(runjags)


#### 2.3 simple model ####
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
# update(jags,2000)
ophtal.sim <- coda.samples(model = jags,
                           variable.names = parameters,
                           n.iter = 5000, 
                           thin = 1)

# Convert osteo.sim into mcmc.list for processing with CODA package
ophtal.mcmc <- as.mcmc.list(ophtal.sim)

# Produce general summary of obtained MCMC sampling
plot(ophtal.mcmc)
summary(ophtal.mcmc)
ophtal.combined <- combine.mcmc(ophtal.mcmc)
HPDinterval(ophtal.combined)

# Specific output obtained from CODA functions
par(mfrow=c(1,3))
traceplot(ophtal.mcmc)
densplot(ophtal.mcmc)

par(mfrow=c(1,1))
acfplot(ophtal.mcmc)
# par(mfrow=c(2,3))
# autocorr.plot(ophtal.mcmc)

par(mfrow=c(1,1))
crosscorr.plot(ophtal.mcmc)

# Convergence tests
gelman.diag(ophtal.mcmc)
gelman.plot(ophtal.mcmc,ask=FALSE)


#### 2.3 Simple Model with prediction ####
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
    
    lutein_2003_new ~ dnorm(mu_new, tau)
    mu_new <- beta0+beta1*5
  }", file="ophtalmology_model_pred.txt")


# Parameters to monitor
parameters <- c("beta0", "beta1", "sigma2", "lutein_2003_new")

# Running JAGS:
jags <- jags.model(file ="ophtalmology_model_pred.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 2)
# update(jags,2000)
ophtal.sim <- coda.samples(model = jags,
                           variable.names = parameters,
                           n.iter = 5000, 
                           thin = 1)

# Convert osteo.sim into mcmc.list for processing with CODA package
ophtal.mcmc_pred <- as.mcmc.list(ophtal.sim)

# Produce general summary of obtained MCMC sampling
plot(ophtal.mcmc_pred)
summary(ophtal.mcmc_pred)
ophtal.combined <- combine.mcmc(ophtal.mcmc_pred)
HPDinterval(ophtal.combined)

# Specific output obtained from CODA functions
par(mfrow=c(1,4))
traceplot(ophtal.mcmc_pred)
densplot(ophtal.mcmc_pred)

par(mfrow=c(1,1))
acfplot(ophtal.mcmc_pred)
# par(mfrow=c(2,3))
# autocorr.plot(ophtal.mcmc)

# Convergence tests
gelman.diag(ophtal.mcmc_pred)
gelman.plot(ophtal.mcmc_pred,ask=FALSE)


#### 2.3 Longer model with thinning, burnin ####
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
plot(ophtal.mcmc_semicomp)
summary(ophtal.mcmc_semicomp)
ophtal.combined <- combine.mcmc(ophtal.mcmc_semicomp)
HPDinterval(ophtal.combined)

# Specific output obtained from CODA functions
par(mfrow=c(1,3))
traceplot(ophtal.mcmc_semicomp)
densplot(ophtal.mcmc_semicomp)

par(mfrow=c(1,1))
acfplot(ophtal.mcmc_semicomp)
# par(mfrow=c(2,3))
# autocorr.plot(ophtal.mcmc_semicomp)

par(mfrow=c(1,1))
crosscorr.plot(ophtal.mcmc_semicomp)

# Convergence tests
gelman.diag(ophtal.mcmc_semicomp)
gelman.plot(ophtal.mcmc_semicomp,ask=FALSE)


#### 2.3 Longer model with thinning, burnin with prediction ####
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
    
    lutein_2003_new ~ dnorm(mu_new, tau)
    mu_new <- beta0+beta1*5
  }", file="ophtalmology_model_semicomp_pred.txt")


# Parameters to monitor
parameters <- c("beta0", "beta1", "sigma2", "lutein_2003_new")

# Running JAGS:
jags <- jags.model(file ="ophtalmology_model_semicomp_pred.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 2)
update(jags,2000)
ophtal.sim <- coda.samples(model = jags,
                           variable.names = parameters,
                           n.iter = 90000, 
                           thin = 10)

# Convert osteo.sim into mcmc.list for processing with CODA package
ophtal.mcmc_semicomp_pred <- as.mcmc.list(ophtal.sim)

# Produce general summary of obtained MCMC sampling
plot(ophtal.mcmc_semicomp_pred)
summary(ophtal.mcmc_semicomp_pred)
ophtal.combined <- combine.mcmc(ophtal.mcmc_semicomp_pred)
HPDinterval(ophtal.combined)

# Specific output obtained from CODA functions
par(mfrow=c(1,4))
traceplot(ophtal.mcmc_semicomp_pred)
densplot(ophtal.mcmc_semicomp_pred)

par(mfrow=c(1,1))
acfplot(ophtal.mcmc_semicomp_pred)
# par(mfrow=c(2,3))
# autocorr.plot(ophtal.mcmc_semicomp_pred)

par(mfrow=c(1,1))
crosscorr.plot(ophtal.mcmc_semicomp_pred)

# Convergence tests
gelman.diag(ophtal.mcmc_semicomp_pred)
gelman.plot(ophtal.mcmc_semicomp_pred,ask=FALSE)


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
plot(ophtal.mcmc_complex)
summary(ophtal.mcmc_complex)
ophtal.combined <- combine.mcmc(ophtal.mcmc_complex)
HPDinterval(ophtal.combined)

# Specific output obtained from CODA functions
par(mfrow=c(1,3))
traceplot(ophtal.mcmc_complex)
densplot(ophtal.mcmc_complex)

par(mfrow=c(1,1))
acfplot(ophtal.mcmc_complex)
# par(mfrow=c(2,3))
# autocorr.plot(ophtal.mcmc_complex)

par(mfrow=c(1,1))
crosscorr.plot(ophtal.mcmc_complex)

# Convergence tests
gelman.diag(ophtal.mcmc_complex)
gelman.plot(ophtal.mcmc_complex,ask=FALSE)


#### 2.3 Complex Model with prediction ####
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
parameters <- c("beta0", "beta1", "sigma2", "lutein_2003_new")

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
plot(ophtal.mcmc_complex_pred)
summary(ophtal.mcmc_complex_pred)
ophtal.combined <- combine.mcmc(ophtal.mcmc_complex_pred)
HPDinterval(ophtal.combined)

# Specific output obtained from CODA functions
par(mfrow=c(1,4))
traceplot(ophtal.mcmc_complex_pred)
densplot(ophtal.mcmc_complex_pred)

par(mfrow=c(1,1))
acfplot(ophtal.mcmc_complex_pred)
# par(mfrow=c(2,3))
# autocorr.plot(ophtal.mcmc_complex_pred)

# Convergence tests
gelman.diag(ophtal.mcmc_complex_pred)
gelman.plot(ophtal.mcmc_complex_pred,ask=FALSE)


#### Comparison of model results (should be around the same) ####
summary(ophtal.mcmc)
summary(ophtal.mcmc_semicomp)
summary(ophtal.mcmc_complex)
summary(ophtal.mcmc_pred)
summary(ophtal.mcmc_semicomp_pred)
summary(ophtal.mcmc_complex_pred)