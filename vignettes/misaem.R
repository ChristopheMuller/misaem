## -----------------------------------------------------------------------------
library(misaem)

## -----------------------------------------------------------------------------
ampute <- function(X.full, prc.NA){

  n <- dim(X.full)[1]
  d <- dim(X.full)[2]

  # Create Mask Matrix
  M <- matrix(0, nrow = n, ncol = d)
  for (j in 1:d) {
    M[, j] <- rbinom(n, 1, prc.NA)
  }
  while (any(rowSums(M) == d)) {
    M.temp <- M[rowSums(M) < d, ]
    left.to.add <- n - nrow(M.temp)
    M.add <- matrix(0, nrow = left.to.add, ncol = d)
    for (j in 1:d) {
      M.add[, j] <- rbinom(left.to.add, 1, prc.NA)
    }
    M <- rbind(M.temp, M.add)
  }


  X.obs <- X.full
  X.obs[M == 1] <- NA

  return(X.obs)

}


## -----------------------------------------------------------------------------
set.seed(1)
n <- 50 # number of rows
p <- 2 # number of explanatory variables

# Generate complete design matrix
library(MASS)
mu.X <- c(1, 1) 
Sigma.X <- matrix(c(1, 1, 1, 4), nrow = 2)
X.complete <- mvrnorm(n, mu.X, Sigma.X)

# Generate response
b <- c(2, 3, -1) # regression coefficient
sigma.eps <- 0.25 # noise variance
y <- cbind(rep(1, n), X.complete) %*% b + rnorm(n, 0, sigma.eps)

## -----------------------------------------------------------------------------
# Add missing values
X.obs <- ampute(X.complete, 0.15)
y.obs <- y

## -----------------------------------------------------------------------------
head(X.obs)

## -----------------------------------------------------------------------------
sum(is.na(X.obs))/(n*p) 

## -----------------------------------------------------------------------------
# Estimate regression using EM with NA 
df.obs = data.frame(y, X.obs)
miss.list = miss.lm(y~., data = df.obs)

## -----------------------------------------------------------------------------
print(miss.list)
print(summary(miss.list)) 
summary(miss.list)$coef 

## -----------------------------------------------------------------------------
# Estimate regression using self-defined parameters
miss.list2 = miss.lm(y~., data = df.obs, print_iter = FALSE, maxruns = 500, tol_em = 1e-4)
print(miss.list2)

## -----------------------------------------------------------------------------
# Add null variable with NA
X.null <- mvrnorm(n, 1, 1)
patterns <- runif(n)<0.15 # missing completely at random
X.null[patterns] <- NA
X.obs.null <- cbind.data.frame(X.obs, X.null)

# Without model selection
df.obs.null = data.frame(y, X.obs.null)
miss.list.null = miss.lm(y~., data = df.obs.null)
print(miss.list.null)

# Model selection
miss.model = miss.lm.model.select(y, X.obs.null)
print(miss.model)

## -----------------------------------------------------------------------------
# Prediction
# Generate dataset
set.seed(200)
nt <- 20  # number of new observations
Xt <- mvrnorm(nt, mu.X, Sigma.X)

# Add missing values
Xt.obs <- ampute(Xt, 0.15)

## -----------------------------------------------------------------------------
#train with NA + test no NA
miss.comptest.pred = predict(miss.list2, data.frame(Xt), seed = 100)
print(miss.comptest.pred)

## -----------------------------------------------------------------------------
#both train & test with NA
miss.pred = predict(miss.list2, data.frame(Xt.obs), seed = 100)
print(miss.pred)

## -----------------------------------------------------------------------------
# Generate dataset
set.seed(200)
n <- 500  # number of subjects
p <- 5     # number of explanatory variables
mu.star <- 1:p  #rep(0,p)  # mean of the explanatory variables
sd <- 1:p # rep(1,p) # standard deviations
C <- matrix(c(   # correlation matrix
1,   0.8, 0,   0,   0,
0.8, 1,   0,   0,   0,
0,   0,   1,   0.3, 0.6,
0,   0,   0.3, 1,   0.7,
0,   0,   0.6, 0.7, 1), nrow=p)
Sigma.star <- diag(sd)%*%C%*%diag(sd) # covariance matrix
beta.star <- c(1, -1, 1, 1, -1) # coefficients
beta0.star <- 0  # intercept
beta.true = c(beta0.star,beta.star)

# Design matrix
X.complete <- matrix(rnorm(n*p), nrow=n)%*%chol(Sigma.star)+
              matrix(rep(mu.star,n), nrow=n, byrow = TRUE)

# Reponse vector
p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
y <- as.numeric(runif(n)<p1)

## -----------------------------------------------------------------------------
# Generate missingness
set.seed(200)
p.miss <- 0.10
X.obs <- ampute(X.complete, p.miss)

## -----------------------------------------------------------------------------
head(X.obs)

## -----------------------------------------------------------------------------
df.obs = data.frame(y, X.obs)

#logistic regression with NA 
miss.list = miss.glm(y~., data = df.obs, seed = 100)
print(miss.list)
print(summary(miss.list)) 
summary(miss.list)$coef 

## -----------------------------------------------------------------------------
# Add null variable with NA
X.null <- mvrnorm(n, 1, 1)
patterns <- runif(n)<0.10 # missing completely at random
X.null[patterns] <- NA
X.obs.null <- cbind.data.frame(X.obs, X.null)

# Without model selection
df.obs.null = data.frame(y, X.obs.null)
miss.list.null = miss.glm(y~., data = df.obs.null)
print(miss.list.null)

# model selection for SAEM
miss.model = miss.glm.model.select(y, X.obs.null)
print(miss.model)

## -----------------------------------------------------------------------------
# Generate test set with missingness
set.seed(200)
nt = 100
X.test <- matrix(rnorm(nt*p), nrow=nt)%*%chol(Sigma.star)+
          matrix(rep(mu.star,nt), nrow = nt, byrow = TRUE)

# Generate the test set
p1 <- 1/(1+exp(-X.test%*%beta.star-beta0.star))
y.test <- as.numeric(runif(nt)<p1)

# Generate missingness on test set
p.miss <- 0.10
X.test[runif(nt*p)<p.miss] <- NA

# Prediction on test set
pr.saem <- predict(miss.list, data.frame(X.test))

# Confusion matrix
pred.saem = (pr.saem>0.5)*1
table(y.test,pred.saem )

