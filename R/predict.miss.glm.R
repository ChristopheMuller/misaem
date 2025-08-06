# Predict Method for miss.glm Fits
#'
#' Prediction on test with missing values for the logistic regression model.
#' @param object a fitted object of class inheriting from "miss.glm".
#' @param newdata a data frame in which to look for variables with which to predict. It can contain missing values.
#' @param seed  An integer as a seed set for the random generator.
#' @param method The name of method to deal with missing values in test set. It can be 'map'(maximum a posteriori) or 'impute' (imputation by conditional expectation). Default is 'map'.
#' @param mc.size The number of Monte Carlo samples to use when method is 'map'. Default is 100.
#' @param ... Further arguments passed to or from other methods.
#' @return
#' \item{pr.saem}{The prediction result for logistic regression: the probability of response y=1.}
#' @import mvtnorm stats MASS abind
#' @importFrom methods is
#' @examples
#' # Generate dataset
#' N <- 100  # number of subjects
#' p <- 3     # number of explanatory variables
#' mu.star <- rep(0,p)  # mean of the explanatory variables
#' Sigma.star <- diag(rep(1,p)) # covariance
#' beta.star <- c(1, 1,  0) # coefficients
#' beta0.star <- 0 # intercept
#' beta.true = c(beta0.star,beta.star)
#' X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) +
#'               matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
#' p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
#' y <- as.numeric(runif(N)<p1)
#'
#' # Generate missingness
#' p.miss <- 0.10
#' patterns <- runif(N*p)<p.miss #missing completely at random
#' X.obs <- X.complete
#' X.obs[patterns] <- NA
#'
#' df.obs = data.frame(y,X.obs)
#'
#' # SAEM
#' miss.list = miss.glm(y~., data=df.obs, print_iter=FALSE,seed=100)
#'
#' # Generate new dataset for prediction
#' Nt <- 20
#' Xt <- matrix(rnorm(Nt*p), nrow=Nt)%*%chol(Sigma.star)+
#'   matrix(rep(mu.star,Nt), nrow=Nt, byrow = TRUE)

#' # Generate missingness in new dataset
#' patterns <- runif(Nt*p)<p.miss
#' Xt.obs <- Xt
#' Xt.obs[patterns] <- NA
#'
#' # Prediction with missing values
#' miss.prob = predict(miss.list, data.frame(Xt.obs), method='map')
#' print(miss.prob)
#' @export

predict.miss.glm <- function(object, newdata = NULL, seed = NA, method='map', mc.size=100, ...)
{
  if (!is.na(seed))
    set.seed(seed)

  X.test = newdata
  mu.saem = object$mu.X
  sig2.saem = object$Sig.X
  beta.saem = object$coef


  #judge
  if (is(X.test, "data.frame")){
    X.test <- as.matrix(X.test)
  }
  if (method == "MAP" | method == "Map") {
    method <- "map"
  }
  if (method == "Impute" | method == "IMPUTE") {
    method <- "impute"
  }

  rindic = as.matrix(is.na(X.test))

  if(sum(rindic)!=0){
    if(method=='impute'){
      for(i in 1:dim(X.test)[1]){
        if(sum(rindic[i,])!=0){
          miss_col = which(rindic[i,]==TRUE)
          x2 = X.test[i,-miss_col]
          mu1 = mu.saem[miss_col]
          mu2 = mu.saem[-miss_col]
          sigma11 = sig2.saem[miss_col,miss_col]
          sigma12 = sig2.saem[miss_col,-miss_col]
          sigma22 = sig2.saem[-miss_col,-miss_col]
          sigma21 = sig2.saem[-miss_col,miss_col]
          mu_cond = mu1+sigma12 %*% solve(sigma22) %*% (x2-mu2)
          X.test[i,miss_col] =mu_cond
        }
      }
      tmp <- as.matrix(cbind.data.frame(rep(1,dim(X.test)[1]),X.test)) %*% as.matrix(beta.saem)
      pr.saem <- 1/(1+(1/exp(tmp)))

    }else if(method=='map'){

      patterns = unique(rindic)
      pr.saem = rep(0, nrow(X.test))

      for (i in 1:nrow(patterns)) {
        pattern = patterns[i, ]
        rows_with_pattern = which(apply(rindic, 1, function(x) all(x == pattern)))
        n_pattern = length(rows_with_pattern)
        
        if (sum(pattern) == 0) {
          x_subset = X.test[rows_with_pattern, , drop = FALSE]
          x_with_intercept = cbind(1, x_subset)
          probs = log_reg(y = 1, x = x_with_intercept, beta.saem, iflog = FALSE)
          pr.saem[rows_with_pattern] = probs
        } else {
          miss_col = which(pattern)
          obs_col = which(!pattern)
          n_missing = length(miss_col)
        
          mu1 = mu.saem[miss_col]
          mu2 = mu.saem[obs_col]
          
          sigma11 = sig2.saem[miss_col, miss_col, drop = FALSE]
          sigma12 = sig2.saem[miss_col, obs_col, drop = FALSE]
          sigma22 = sig2.saem[obs_col, obs_col, drop = FALSE]
          sigma21 = sig2.saem[obs_col, miss_col, drop = FALSE]

          inv_sigma22 = solve(sigma22)
          inv_sigma22_sigma21 = inv_sigma22 %*% sigma21

          sigma_cond = sigma11 - sigma12 %*% inv_sigma22_sigma21
          mu_cond_diff = sigma12 %*% inv_sigma22 %*% t(x2 - mu2)
          mu_cond = mu1 + mu_cond_diff
          mu_cond = t(mu_cond)

          x1_samples_list = apply(mu_cond, 1, function(row_mu) {
                rmvnorm(n = mc.size, mean = row_mu, sigma = sigma_cond)
              }, simplify=FALSE)      # dimension: (n_missing x mc.size, n_pattern)
          x1_all = simplify2array(x1_samples_list)

          x_observed_replicated = replicate(mc.size, X.test[rows_with_pattern, , drop = FALSE], simplify = "array") # dim = [n_pattern, ncol(X.test), mc.size]
          for (j in 1:n_missing) {
            x_observed_replicated[, miss_col[j], ] = x1_all[, j, ]
          }
          x_imputed_flat = aperm(x_observed_replicated, c(2, 1, 3)) # dim = [ncol(X.test), n_pattern, mc.size]
          dim(x_imputed_flat) = c(mc.size * n_pattern, dim(x_observed_replicated)[2])
        
          x_imputed_with_intercept = cbind(1, x_imputed_flat)
          
          linear_pred = log_reg(y = 1, x = x_imputed_with_intercept, beta.saem, iflog = FALSE)
          probs_matrix = matrix(linear_pred, nrow = mc.size, ncol = n_pattern)
    
          pr.saem[rows_with_pattern] = colMeans(probs_matrix)
        
        }
      }

      pr.saem = as.matrix(pr.saem)

    } else {
      stop("Error: There is no such method. Method should be 'map' or 'impute'. ")
    }
  } else {
    tmp <- as.matrix(cbind.data.frame(rep(1,dim(X.test)[1]),X.test)) %*% as.matrix(beta.saem)
    pr.saem <- 1/(1+(1/exp(tmp)))
  }


  return(pr.saem)
}
