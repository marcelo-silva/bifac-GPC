##################################################
## bifac-GPC model with logit link function     ##
##################################################

##################################################
# Source code for Application of the manuscript: #
# "A bifactor generalized partial credit model   #
#  with flexible link functions: a novel         #
#  approach in survey data", by                  #
#  Marcelo A. da Silva, Anne C. Huggins-Manley,  #
#  José A. Mazzon, and Jorge L. Bazán            #
##################################################

##### Set the path of files #####
setwd("C:/users/application_JAS")

##### Load R packages #####
library(rstan)    # R package to interact with Stan
library(loo)      # R package to calculate the LOO criterion

##### Functions #####
# Function to calculate the WAIC criterion
waic <- function(stanfit){ # Function to calculate the WAIC criterion
  log_lik <- extract (stanfit, "log_lik")$log_lik
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
  c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
  matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
  colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
  p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
  pointwise=pointwise, total=total, se=se))
}

# Function to calculate the DIC criterion
functionDhat <- function(y, prob) {
  categ <- array(dim=length(Y))
  lenY <- length(Y)
  for (l in 1:lenY) {
    categ[l] <- log(prob[l, y[l]])
  }
  return(-2*sum(categ))
}

# Function to transform multidimensional vector into one-dimensional vector
vectorization <- function(M) {
  dimension <- dim(M)
  nNA <- sum(length(which(is.na(M))))
  index <- 1:(prod(dimension)-nNA)
  if (is.element(NA, M)) {
    idxNA <- which(is.na(M))
    for (i in 1:nNA) 
      index <- append(index, 0, after=(idxNA[i]-1))
  }
  idx <- array(index, dim=dimension)
  Mvec <- as.vector(M)
  MvecFinal <- Mvec[!is.na(Mvec)]
  vector <- list("vec"=MvecFinal, "idx"=idx)
  return(vector)
}

# Function for constructing indices of discrimination parameters
a_indexes <- function(J, D) {
  maxJ <- max(J)
  indexes <- c(1:J[1], rep(0, maxJ-J[1]))
  for (d in 2:D)
    indexes <- c(indexes, (max(indexes)+1):(max(indexes)+J[d]), rep(0, maxJ-J[d]))
  a_idx <- array(indexes, dim=c(maxJ, D))
  return(a_idx)
}

# Function for constructing indices of difficulty parameters
b_indexes <- function(J, D, m) {
  maxJ <- max(J)
  maxm <- max(m)
  for (k in 1:(maxm-1)) {
    for (d in 1:D) {
      if ((k == 1) && (d == 1)) {
        indexes <- c(1:J[d], rep(0, maxJ-J[d]))
      } else {
        indexes <- c(indexes, (max(indexes)+1):(max(indexes)+J[d]), rep(0, maxJ-J[d]))
      }
    }
  }
  b_idx <- array(indexes, dim=c(maxJ, D, maxm-1))
  return(b_idx)
}


##### Settings for the application #####
# Number of individuals
n <- 333

# Number of items in each dimension
J <- c(4,5,4,4,5,4,4)

# Maximum number of items in dimensions
maxJ <- max(J)

# Total number of items
sumJ <- sum(J)

# Number of dimensions
D <- length(J)

# Vector with number of response categories
m <- array(c(6,6,6,6,0,
             6,6,6,6,6,
             6,6,6,6,0,
             6,6,6,6,0,
             6,6,6,6,6,
             6,6,6,6,0,
             6,6,6,6,0), dim=c(maxJ, D))

# Maximum number of response categories
maxm <- max(m)

# Total number of response categories
summ <- sum(m)

##### Mobile banking data #####
# Reading the CSV file with the data
data <- read.csv("dados1.csv", sep=";", header=TRUE)

# Organization of data in a 3-dimensional vector
data3d <- array(data=0, dim=c(n, maxJ, D))
data3d[,,1] <- as.matrix(cbind(data[,1:4], rep(NA, n)))    # Compatibility
data3d[,,2] <- as.matrix(data[,5:9])                       # Relative advantage
data3d[,,3] <- as.matrix(cbind(data[,10:13], rep(NA, n)))  # Visibility
data3d[,,4] <- as.matrix(cbind(data[,14:17], rep(NA, n)))  # Results demonstrability
data3d[,,5] <- as.matrix(data[,18:22])                     # Image
data3d[,,6] <- as.matrix(cbind(data[,23:26], rep(NA, n)))  # Triability
data3d[,,7] <- as.matrix(cbind(data[,27:30], rep(NA, n)))  # Perceived ease of use

##### Vectoring the 3-dimensional vector data3d #####
# Transform 3-dimensional vector into one-dimensional vector
Y_array <- vectorization(data3d)

# Saves the one-dimensional vector
Y <- Y_array$vec

# Save the vector of indices
Y_idx <- Y_array$idx

# Length of the one-dimensional vector
lenY <- length(Y)

##### Constructing the index of the parameters #####
# Discrimination parameters (global and specific)
a_idx <- a_indexes(J, D)

# Difficulty  parameters
b_idx <- b_indexes(J, D, m)

# Global latent trait
theta0_idx <- 1:n

# Specific latent traits
theta_idx <- array(1:(D*n), dim=c(D, n))

##### Settings for executing Stan code (estimates and criteria) #####
# Parameters to be returned
param <- c("a0", "a", "b", "theta0", "theta", "prob", "dev", "log_lik")

#Number of iterations
niter <- 10000

# Number of chains
nchains <- 1

# Stan Function
stanfit <- stan(file   = 'bifacGPClogit.stan',
                pars   = param,
                chains = nchains,
                seed   = 7,
                thin   = 1,
                iter   = niter)

##### Save results in CSV file #####
# Summary of results
results <- summary(stanfit, probs = NA)$summary

# Selecting the model parameters
results_stan <- cbind(results[1:(maxJ*D*(maxm+1)+(D+1)*n), -4])

# Saves estimates to CSV file
write.table(results_stan, file="stanfitLogit.csv", sep=";",
            dec=".", col.names=TRUE, row.names=TRUE)

##### Model comparison criteria #####
# This value is set to -2* logarithm of the likelihood
Dbar <- mean(extract(stanfit, "dev", inc_warmup=FALSE, permuted=FALSE))

# Probabilities according to the bifac-GPC model
prob <- array(colMeans(extract(stanfit, "prob", permuted=FALSE)), 
              dim=c(lenY, max(m)))

# DIC calculation
Dhat <- functionDhat(Y, prob)
pD <- Dbar - Dhat
DIC <- Dbar + pD

# WAIC calculation
WAIC <- waic(stanfit)

# Logarithm of likelihood
log_lik <- extract_log_lik(stanfit)

# LOO calculation
loo <- loo(log_lik)

##### Saves model comparison criteria values in CSV file #####
# Saves criteria values to CSV file
crit_values <- matrix(c(DIC,WAIC$waic,loo$looic),nrow=1)
write.table(crit_values, file="criteriaLogit.csv", sep=";", 
            col.names=c("DIC","WAIC","loo"), 
            dec=".", row.names=FALSE)


##### Settings for executing Stan code (PPMC)#####
# Parameter to be returned
param <- c("y_rep")

# Number of iterations
niter <- 10000

# Number of chains
nchains <- 1

# Stan Function
stanfit <- stan(file   = 'bifacGPClogitrep.stan',
                pars   = param,
                chains = nchains,
                seed   = 7,
                thin   = 1,
                iter   = niter)

##### Calculate Bayesian p-value #####
# Saves the replicates of the response vector y
y_rep <- extract(stanfit, "y_rep", inc_warmup=FALSE, permuted=FALSE)

# Sets the number of replicas to be used (excludes replicas of burn-in)
T <- (niter/2)*nchains

# Declares the frequency vector
F <- array(NA, dim=c(T, 6))

# Saves the frequency in each replica
for (t in 1:T) 
  F[t,] <- table(array(y_rep[t,1,]))

# Calculates the expected frequency
Expectation <- colMeans(F)

# Calculates the observed discrepancy measure
D_obs <- sum((table(Y) - Expectation)^2/Expectation)

# Auxiliary Variable
indicator <- 0

# Declares the discrepancy measure of replicas
D_rep <- array(NA, dim=T)

# Counts the number of D_rep are greater than or equal to D_obs
for (t in 1:T) {
  D_rep[t] <- sum((F[t,] - Expectation)^2/Expectation)
  if (D_rep[t] >= D_obs) {
    indicator <- indicator + 1
  }
}

# Calculates the Bayesian p-value by the ratio "indicator"
pb <- indicator/T

##### Saves Bayesian p-value in CSV file #####
write.table(pb, file="pvalueLogit.csv", sep=";", 
            col.names=c("Bayesian p-value"), dec=".")
