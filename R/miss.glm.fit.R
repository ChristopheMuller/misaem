#' Fitting Logistic Regression Models with Missing Values
#'
#' This function is used inside \code{miss.glm} to fit logistic regression model with missing values, by algorithm SAEM.
#' @param x design matrix with missingness \eqn{N \times p}{N * p}.
#' @param y response vector \eqn{N \times 1}{N * 1}.
#' @param control a list of parameters for controlling the fitting process. For \code{miss.glm.fit} this is passed to \code{\link{miss.glm.control}}.
#' @return a list with following components:
#' \item{coefficients}{Estimated \eqn{\beta}{\beta}.}
#' \item{ll}{Observed log-likelihood.}
#' \item{var.covar}{Variance-covariance matrix for estimated parameters.}
#' \item{s.err}{Standard error for estimated parameters.}
#' \item{mu.X}{Estimated \eqn{\mu}{\mu}.}
#' \item{Sig.X}{Estimated \eqn{\Sigma}{\Sigma}.}
#' @import mvtnorm stats glmnet
#' @examples
#' ## For examples see example(miss.glm)

miss.glm.fit <- function (x, y, control = list()) {
  control <- do.call("miss.glm.control", control)
  xnames <- dimnames(x)[[2L]]
  x <- as.matrix(x[,-1])

  maxruns <- control$maxruns
  tol_em <- control$tol_em
  nmcmc <- control$nmcmc
  tau <- control$tau
  k1 <- control$k1
  seed <- control$seed
  print_iter <- control$print_iter
  var_cal <- control$var_cal
  ll_obs_cal <- control$ll_obs_cal
  subsets <- control$subsets
  save_trace <- control$save_trace
  alphaReg <- control$alpha
  lambda <- control$lambda

  if (!is.na(seed))
    set.seed(seed)

  p <- ncol(x)

  if (is.na(subsets[1]))
    subsets <- 1:p

  if (sum(subsets %in% 1:p) < length(subsets)) {
    stop("Error: index of selected variables must be in the range of covariates.")
  }
  if (length(unique(subsets)) != length(subsets)) {
    stop("Error: index of selected variables must not be repeated.")
  }

  if (any(rowSums(is.na(x)) == p)) {
    i_allNA <- which(rowSums(is.na(x)) == p)
    x <- x[-i_allNA, ]
    y <- y[-i_allNA]
  }
  if (any(is.na(y))) {
    i_YNA <- which(is.na(y))
    x <- x[-i_YNA, ]
    y <- y[-i_YNA]
  }

  if (sum(sapply(x, is.numeric)) < ncol(x)) {
    stop("Error: the variables should be numeric.")
  }
  if (sum(y == 1) + sum(y == 0) < nrow(x)) {
    stop("Error: y must be coded by 0 or 1, and there is no missing data in y.")
  }

  n <- length(y)

  rindic <- as.matrix(is.na(x))
  if(sum(rindic)>0){
    whichcolmissing <- (1:ncol(rindic))[apply(rindic,2,sum)>0]
    missingcols <- length(whichcolmissing)
  }
  if(sum(rindic)==0){
    whichcolmissing <- integer(0)
    missingcols <- 0
  }

  missingcols <- as.integer(sum(colSums(rindic) > 0))

  if (missingcols > 0) {
    k <- 0
    cstop <- 0.1
    seqbeta <- matrix(NA, nrow = ncol(x) + 1, ncol = (maxruns + 1))
    seqbeta_avg <- matrix(NA, nrow = ncol(x) + 1, ncol = (maxruns + 1))

    X.sim <- X.mean <- initialize_df(x = x, method = control$init_method, seed = seed)

    mu <- colMeans(X.mean)
    Sigma <- var(X.mean) * (n - 1) / n
    beta <- rep(0, p + 1)
    if (lambda == 0){
      beta[c(1,subsets+1)]= glm(y~ X.mean[,subsets],family=binomial(link='logit'))$coef
    } else {
      g <- glmnet(X.mean[,subsets], y, family = "binomial", alpha = alphaReg, lambda = lambda)
      beta[subsets+1] <- g$beta[subsets, 1]
      beta[1] <- g$a0[1]
    }

    if(print_iter==TRUE){
      cat(sprintf('Iteration of SAEM: \n'))
    }

    pattern_ids <- apply(rindic, 1, paste, collapse = "_")
    grouped_rows <- split(seq_len(n), pattern_ids)

    while ((cstop > tol_em && k < maxruns) || (k < 20)) {
      k <- k + 1
      beta.old <- beta

      gamma <- if (k < k1) 1 else 1 / (k - (k1 - 1))^tau

      S.inv <- chol2inv(chol(Sigma))

      for (rows_with_pattern in grouped_rows) {
        patt <- rindic[rows_with_pattern[1], ]
        jna <- which(patt)
        njna <- length(jna)

        if (njna == 0) {
          next
        }

        n_pattern <- length(rows_with_pattern)
        jobs <- which(!patt)

        S.inv_MM <- S.inv[jna, jna, drop = FALSE]
        Oi <- chol2inv(chol(S.inv_MM))
        chol_Oi <- chol(Oi)

        if (njna < p) {
          X_obs <- X.sim[rows_with_pattern, jobs, drop = FALSE]
          S.inv_MO <- S.inv[jna, jobs, drop = FALSE]
          delta_X <- t(X_obs) - mu[jobs]
          adjustment <- t(Oi %*% S.inv_MO %*% delta_X)
          mi <- matrix(mu[jna], nrow = n_pattern, ncol = njna, byrow = TRUE) - adjustment
          lobs <- beta[1] + X_obs %*% beta[jobs + 1]
        } else {
          mi <- matrix(mu[jna], nrow = n_pattern, ncol = njna, byrow = TRUE)
          lobs <- rep(beta[1], n_pattern)
        }

        xina <- X.sim[rows_with_pattern, jna, drop = FALSE]
        betana <- beta[jna + 1]
        y_pattern <- y[rows_with_pattern]
        is_y1 <- (y_pattern == 1)

        for (m in 1:nmcmc) {
          rand_norm <- matrix(rnorm(n_pattern * njna), nrow = n_pattern, ncol = njna)
          xina.c <- mi + rand_norm %*% chol_Oi

          current_logit <- xina %*% betana + lobs
          candidate_logit <- xina.c %*% betana + lobs

          num_y1 <- 1 + exp(-current_logit)
          den_y1 <- 1 + exp(-candidate_logit)
          ratio_y1 <- num_y1 / den_y1
          fix_idx_y1 <- is.nan(ratio_y1)
          if (any(fix_idx_y1)) {
            delta <- candidate_logit[fix_idx_y1] - current_logit[fix_idx_y1]
            new_ratio <- exp(delta)
            new_ratio[new_ratio == Inf] <- 1.01
            ratio_y1[fix_idx_y1] <- new_ratio
          }

          num_y0 <- 1 + exp(current_logit)
          den_y0 <- 1 + exp(candidate_logit)
          ratio_y0 <- num_y0 / den_y0
          fix_idx_y0 <- is.nan(ratio_y0)
          if (any(fix_idx_y0)) {
            delta <- current_logit[fix_idx_y0] - candidate_logit[fix_idx_y0]
            new_ratio <- exp(delta)
            new_ratio[new_ratio == Inf] <- 1.01
            ratio_y0[fix_idx_y0] <- new_ratio
          }

          alpha <- ifelse(is_y1, ratio_y1, ratio_y0)

          accepted <- runif(n_pattern) < alpha
          if (any(accepted)) {
            xina[accepted, ] <- xina.c[accepted, , drop = FALSE]
          }
        }
        X.sim[rows_with_pattern, jna] <- xina
      }
      beta_new <- rep(0,p+1)
      if (lambda == 0){
        beta_new[c(1,subsets+1)] <- glm(y~ X.sim[,subsets],family=binomial(link='logit'))$coef
      } else {
        g <- glmnet(X.sim[,subsets], y, family = "binomial", alpha = alphaReg, lambda = lambda)
        beta_new[subsets+1] <- g$beta[subsets, 1]
        beta_new[1] <- g$a0[1]
      }
      
      beta <- (1-gamma)*beta + gamma*beta_new
      cstop <- sum((beta-beta.old)^2)

      mu <- (1 - gamma) * mu + gamma * colMeans(X.sim)
      Sigma <- (1 - gamma) * Sigma + gamma * (cov(X.sim) * (n - 1) / n)

      seqbeta[, k] <- beta_new
      seqbeta_avg[, k] <- beta.old

      if (print_iter && k %% 50 == 0) {
        cat(sprintf('%i... ', k))
      }
    }
    var_obs = ll = std_obs =NULL
    if(var_cal==TRUE){
      var_obs = louis_lr_saem(beta,mu,Sigma,y,x,subsets,rindic,whichcolmissing,mc.size=100)
      std_obs <- sqrt(diag(var_obs))
      if (any(is.nan(std_obs))){
        warning("Instabilities in standard errors.")
      }
    }
    if(ll_obs_cal==TRUE){
      ll = likelihood_saem(beta,mu,Sigma,y,x,rindic,whichcolmissing,mc.size=100)
    }

  } else{ # no missing data

    x <- as.matrix(x)
    
    x_sub <- x[, subsets, drop = FALSE]
    
    data.complete <- data.frame(y=y, x_sub)
    mu <- apply(x,2,mean)
    Sigma <- var(x)*(n-1)/n
    beta <- rep(0, p + 1)
    
    if (lambda == 0){
      model.complete <- glm(y ~ ., family=binomial(link='logit'), data=data.complete)
      beta[c(1, subsets + 1)] <- model.complete$coefficients
    } else {
      model.complete <- glmnet(x_sub, y, family = "binomial", alpha = alphaReg, lambda = lambda)
      beta[subsets+1] <- as.vector(model.complete$beta)
      beta[1] <- model.complete$a0
    }
    
    var_obs = ll = std_obs = seqbeta_avg = seqbeta = NULL
    
    if(var_cal==TRUE){
      P <- predict(model.complete, newx = x_sub, type = "response")
      P <- as.vector(P)

      W <- diag(P*(1-P))
      X_design <- cbind(1, x_sub)
      
      var_obs_sub <- solve(t(X_design) %*% W %*% X_design)
      
      var_obs <- matrix(0, p + 1, p + 1)
      var_obs[c(1, subsets + 1), c(1, subsets + 1)] <- var_obs_sub
      std_obs <- sqrt(diag(var_obs))
    }
    
    if(ll_obs_cal==TRUE){
      ll = likelihood_saem(beta,mu,Sigma,y,x,rindic,whichcolmissing,mc.size=100)
    }
  }

  beta = beta[c(1,subsets+1)]
  names(beta) <- xnames[c(1,subsets+1)]
  names(std_obs) <- xnames[c(1,subsets+1)]
  
  result <- list(coefficients = beta, s.err = std_obs, var.covar = var_obs, ll = ll[1,1], Sig.X = Sigma, mu.X = mu)
  
  if (save_trace && missingcols > 0) {
    result$trace <- seqbeta[c(1,subsets+1), 1:k, drop = FALSE]
    result$trace_avg <- seqbeta_avg[c(1,subsets+1), 1:k, drop = FALSE]
  }
  
  return(result)
}
