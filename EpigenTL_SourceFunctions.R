


################################################################################
############################# Helper Functions #################################


# this function aggregates the coefficients from target and auxiliary data
# using the error in the target data. 
# this is called in the background, but is not being called directly by users.

## B corresponds to the coefficient vector
## X is the matrix of covariates 
## y is the outcome vector 
## N_vector is a vector of number of observations in the target and aux datasets in that order. 
coef.aggr <- function(B, X, y, N_vector){
  
  if(sum(B == 0) == ncol(B)*nrow(B)){ # if all coefficients are zero, just return zero
    return(rep(0,nrow(B)))
  }
  p <- nrow(B)
  K <- ncol(B)
  colnames(B) <- NULL
  
    # Take the difference in target y and predicted y from each aux data coefficients 
    y0hatk <- -log(colSums(y[1:N_vector[1]] - X[1:N_vector[1], ] %*% B)^2)
    theta.hat <- exp(y0hatk)
    theta.hat = theta.hat / sum(theta.hat)
    # weights by fraction of total squared error

    # multiply betak by weights for each kth aux 
    beta <- as.numeric(B%*%theta.hat) 

  list(theta = theta.hat, beta = beta)

}



# this tells the TL functions where the first and last observation is for 
# each dataset (target, aux 1, ..., aux k) based on the sample sizes in n.vec. 
## n.vec is the vector of number of observations in the target and aux datasets in that order. 
## k is index of the dataset we are interested in extracting values from. 
ind.set <- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k==1){
      ind.re<-c(ind.re,1: n.vec[1])
    }else{
      ind.re<- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
  }
  ind.re
}


# Not in from dplyr 
`%!in%` = Negate(`%in%`) 





################################################################################
########################### Background Functions ###############################


## This is the actual Transfer Learning Lasso loop called inside both 
# TL.Lasso.EstA0 and TL.Lasso.Oracles
## X is the matrix of covariates 
## y is the outcome vector 
## A0 is the informative Aux Set. It is either NULL or estimated
## n.vec is a vector of number of observations in the target and aux datasets in that order. 
## lam.const is whether we are calculating the optimal lambda via CV in each informative aux set 
## or the constant value if using the methods outlined in TransLasso paper. 
TL_Lasso <- function(X, y, A0, n.vec, lam.const=NULL, lam.const_1se = NULL, ...){
  
  p <- ncol(X)
  size.A0 <- length(A0) # set to NULL so its 0
  
  if(size.A0 > 0){ # only for Aux data, otherwise SKIP to below
    
    ind.kA <- ind.set(n.vec, c(1, A0+1))
    ind.1 <- 1:n.vec[1] # vector of all values to build initial model.
    
    y.A <- y[ind.kA] 
    
    # if null, CV done for each Informative Set and both min and 1se lambda kept.
    if(is.null(lam.const)){ # this gets run on first run because we have it set to NULL
      # does its own grid search for lambda
      cv.init<-cv.glmnet(X[ind.kA,], y.A, nfolds=8)
      #### change here for lambda.1se
      #### now, it will just take whatever the best value was and calculate the constant. 
      lam.const <- cv.init$lambda.min/sqrt(2*log(p)/length(ind.kA))
      lam.const_1se <- cv.init$lambda.1se/sqrt(2*log(p)/length(ind.kA))
    }
    if(!is.null(lam.const) & is.null(lam.const_1se)){
      lam.const_1se <- lam.const
    }
    
    # w.kA = coefficients from Xk predicts Yk
    w.kA_min <- as.numeric(glmnet(X[ind.kA,], y.A, lambda=lam.const*sqrt(2*log(p)/length(ind.kA)))$beta)
    w.kA_1se <- as.numeric(glmnet(X[ind.kA,], y.A, lambda=lam.const_1se*sqrt(2*log(p)/length(ind.kA)))$beta)
    
    # w.k coefficient thresholding
    w.kA_min_halflambda <- w.kA_min*(abs(w.kA_min)>=0.5*lam.const*sqrt(2*log(p)/length(ind.kA))) 
    w.kA_min_lambda <- w.kA_min*(abs(w.kA_min)>=lam.const*sqrt(2*log(p)/length(ind.kA))) 
    
    w.kA_1se_halflambda <- w.kA_1se*(abs(w.kA_1se)>=0.5*lam.const_1se*sqrt(2*log(p)/length(ind.kA))) 
    w.kA_1se_lambda <- w.kA_1se*(abs(w.kA_1se)>=lam.const_1se*sqrt(2*log(p)/length(ind.kA))) 

    # build model in target where outcome is what is left after taking the yhat from kth model. 
    # delta.kA = coefficients from X0 predicts (Y0 - Y0hat from w.kA)
    delta.kA_min <- as.numeric(glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA_min, 
                                      lambda=lam.const*sqrt(2*log(p)/length(ind.1)))$beta)
    delta.kA_1se <- as.numeric(glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA_1se, lambda=lam.const_1se*sqrt(2*log(p)/length(ind.1)))$beta)

    # delta.k coefficient thresholding
    delta.kA_min_halflambda <- delta.kA_min*(abs(delta.kA_min)>=0.5*lam.const*sqrt(2*log(p)/length(ind.1))) 
    delta.kA_min_lambda <- delta.kA_min*(abs(delta.kA_min)>=lam.const*sqrt(2*log(p)/length(ind.1))) 
    
    delta.kA_1se_halflambda <- delta.kA_1se*(abs(delta.kA_1se)>=0.5*lam.const_1se*sqrt(2*log(p)/length(ind.1))) 
    delta.kA_1se_lambda <- delta.kA_1se*(abs(delta.kA_1se)>=lam.const_1se*sqrt(2*log(p)/length(ind.1))) 

    # final beta coefficients (the kth lasso coefficients are the weights from kth lasso model + 
    # kth lasso on target)
    
    # no thresholding
    beta.kA_min <- w.kA_min + delta.kA_min 
    beta.kA_1se <- w.kA_1se + delta.kA_1se
    
    # half lambda thresholding
    beta.kA_min_halflambda <- w.kA_min_halflambda + delta.kA_min_halflambda
    beta.kA_min_lambda <- w.kA_min_lambda + delta.kA_min_lambda
    
    # lambda thresholding
    beta.kA_1se_halflambda <- w.kA_1se_halflambda + delta.kA_1se_halflambda
    beta.kA_1se_lambda <- w.kA_1se_lambda + delta.kA_1se_lambda
    
    lam.const=NULL # reset lambda, but I think it actually skips the initialization for the 
    # first auxillary dataset because we supply the lambda constant. 
    
    # output all coefficients
    # recall if constant lambda was specified, min and 1se will be identical. 
    list(beta.kA_min = as.numeric(beta.kA_min), w.kA_min=w.kA_min, 
       beta.kA_1se = as.numeric(beta.kA_1se), w.kA_1se=w.kA_1se,
       beta.kA_min_halflambda = as.numeric(beta.kA_min_halflambda), w.kA_min_halflambda=w.kA_min_halflambda, 
       beta.kA_1se_halflambda = as.numeric(beta.kA_1se_halflambda), w.kA_1se_halflambda=w.kA_1se_halflambda,
       beta.kA_min_lambda = as.numeric(beta.kA_min_lambda), w.kA_min_lambda=w.kA_min_lambda, 
       beta.kA_1se_lambda = as.numeric(beta.kA_1se_lambda), w.kA_1se_lambda=w.kA_1se_lambda,
       lam.const=lam.const, lam.const_1se=lam.const_1se)
  }else{ # end of if(size.A0 > 0)
    
    # BEGINNING CODE FOR INITIAL / TARGET DATA
    cv.init <- cv.glmnet(X[1:n.vec[1],], y[1:n.vec[1]], nfolds=8)
    
    # When constant lambda selected, min lambda is used. 
    lam.const <- cv.init$lambda.min/sqrt(2*log(p)/n.vec[1])
    lam.const_1se <- cv.init$lambda.1se/sqrt(2*log(p)/n.vec[1])
    
    # extract coefficients
    beta.kA_min <- predict(cv.init, s='lambda.min', type='coefficients')[-1] # we are excluding the intercept...
    beta.kA_1se <- predict(cv.init, s='lambda.1se', type='coefficients')[-1] # we are excluding the intercept...
    w.kA_min <- w.kA_1se <- NA
    
    list(beta.kA_min = as.numeric(beta.kA_min), w.kA_min=w.kA_min, 
       beta.kA_1se = as.numeric(beta.kA_1se), w.kA_1se=w.kA_1se,
       lam.const=lam.const, lam.const_1se = lam.const_1se)
  }
  
}




## This performs Oracle Transfer Learning Lasso meaning the set of auxiliary 
# datasets is specified to either run all at once (Oracle) or 1 dataset at a time
# (Oracle 1df). 
## X is the matrix of covariates 
## y is the outcome vector 
## n.vec is a vector of number of observations in the target and aux datasets in that order. 
## AuxInformation is either "Oracle" or "Oracle 1df". 
## LambdaType is either "Constant" or "CV" to indicate if the lambda calculated 
## from the target dataset should be used and adjusted based on the aux dataset size ("Constant")
## or whether the optimal lambda should be calculated via cross validation ("CV") for
## each set of auxiliary information. 
TransLasso.Oracles <- function(X, y, n.vec, AuxInformation = "Oracle 1df", 
                                  LambdaType = "CV", ...) {
  
  M = length(n.vec)-1
  
  p <- ncol(X)
  ind.1 <- ind.set(n.vec, 1) # make row indices of where these observations are for target
 
  Tset <- list()
  
  if(AuxInformation == "Oracle"){
      # make Tset actually just all the aux datasets (for oracle all at once)
      # aux datasets start at 2.
      Tset[[1]] <- c(1:M) 
  }
  # take 1 aux dataset at a time, noting that we index by dataset before.
  if(AuxInformation == "Oracle 1df"){
    for(kk in 1:M){ #use Rhat as the selection rule
     Tset[[kk]] <- kk
    } # the sets of aux datasets to take for each ranking of datasets to include. 
  }
  
  k0 = length(Tset)
  Tset <- unique(Tset)

  beta.T_min <- beta.T_min_lambda <- beta.T_min_halflambda <- list()
  beta.T_1se <- beta.T_1se_lambda <- beta.T_1se_halflambda <- list()
  
  # Lasso on Target Data only 
  init.re <- TL_Lasso(X=X, y=y, A0=NULL, n.vec=n.vec)
  
  beta.T_min[[1]] <- init.re$beta.kA_min 
  beta.T_1se[[1]] <- init.re$beta.kA_1se
  
  # if constant lambda specified, it is here. 
  c1_lambda_const <- init.re$lam.const
  c1_lambda_const_1se <- init.re$lam.const_1se
  
  og_lasso_coef_min <- init.re$beta.kA_min
  og_lasso_coef_1se <- init.re$beta.kA_1se
  
  
  beta.T_min_lambda <- beta.T_min_halflambda <- beta.T_min
  beta.T_1se_lambda <- beta.T_1se_halflambda <- beta.T_1se
  
  beta.pool.T_min <- beta.pool.T_min_lambda <- beta.pool.T_min_halflambda <- beta.T_min ## another method for comparison
  beta.pool.T_1se <- beta.pool.T_1se_lambda <- beta.pool.T_1se_halflambda <- beta.T_1se ## another method for comparison
  
  # go through TL Lasso for each informative set 
  for(kk in 1:length(Tset)){
    T.k <- Tset[[kk]]
    
    # changed function call to lam.const = NULL to do Aux CV
    if(LambdaType == "CV"){
      re.k <- TL_Lasso(X=X, y=y, A0=T.k, n.vec=n.vec, lam.const = NULL)
    }else{ # constant lambda specification
      re.k <- TL_Lasso(X=X, y=y, A0=T.k, n.vec=n.vec, lam.const = c1_lambda_const,
                       lam.const_1se = c1_lambda_const_1se)
    }
    
    # extract coefficients for each informative auxiliary set
    beta.T_min[[kk+1]] <- re.k$beta.kA_min 
    beta.pool.T_min[[kk+1]] <- re.k$w.kA_min
    
    beta.T_1se[[kk+1]] <- re.k$beta.kA_1se 
    beta.pool.T_1se[[kk+1]] <- re.k$w.kA_1se
    
    beta.T_min_halflambda[[kk+1]] <- re.k$beta.kA_min_halflambda 
    beta.pool.T_min_halflambda[[kk+1]] <- re.k$w.kA_min_halflambda
    
    beta.T_1se_lambda[[kk+1]] <- re.k$beta.kA_1se_lambda 
    beta.pool.T_1se_lambda[[kk+1]] <- re.k$w.kA_1se_lambda
  }
  
  beta.T_min <- beta.T_min[!duplicated((beta.T_min))] 
  beta.T_min <- as.matrix(as.data.frame(beta.T_min)) 
  
  beta.T_1se <- beta.T_1se[!duplicated((beta.T_1se))] 
  beta.T_1se <- as.matrix(as.data.frame(beta.T_1se)) 
  
  beta.T_min_halflambda <- beta.T_min_halflambda[!duplicated((beta.T_min_halflambda))] 
  beta.T_min_halflambda <- as.matrix(as.data.frame(beta.T_min_halflambda)) 
  beta.T_min_lambda <- beta.T_min_lambda[!duplicated((beta.T_min_lambda))] 
  beta.T_min_lambda <- as.matrix(as.data.frame(beta.T_min_lambda)) 
  
  beta.T_1se_halflambda <- beta.T_1se_halflambda[!duplicated((beta.T_1se_halflambda))] 
  beta.T_1se_halflambda <- as.matrix(as.data.frame(beta.T_1se_halflambda))
  beta.T_1se_lambda <- beta.T_1se_lambda[!duplicated((beta.T_1se_lambda))] 
  beta.T_1se_lambda <- as.matrix(as.data.frame(beta.T_1se_lambda)) 
  
  ## aggregate coefficients using squared error. 
  ## aggregate w.kA and delta.kA
  agg.re1_min <- coef.aggr(B= beta.T_min, X = X, y = y, N_vector = n.vec)
  agg.re1_1se <- coef.aggr(B= beta.T_1se, X = X, y = y, N_vector = n.vec)
  
  agg.re1_min_halflambda <- coef.aggr(B= beta.T_min_halflambda, X = X, y = y, N_vector = n.vec)
  agg.re1_1se_halflambda <- coef.aggr(B= beta.T_1se_halflambda, X = X, y = y, N_vector = n.vec)
  
  agg.re1_min_lambda <- coef.aggr(B= beta.T_min_lambda, X = X, y = y, N_vector = n.vec)
  agg.re1_1se_lambda <- coef.aggr(B= beta.T_1se_lambda, X = X, y = y, N_vector = n.vec)
  
  
  
  return(list(beta.hat_min = agg.re1_min$beta, theta.hat_min = agg.re1_min$theta, 
              beta.hat_1se = agg.re1_1se$beta, theta.hat_1se = agg.re1_1se$theta,
              beta.hat_min_halflambda = agg.re1_min_halflambda$beta, theta.hat_min_halflambda = agg.re1_min_halflambda$theta, 
              beta.hat_1se_halflambda = agg.re1_1se_halflambda$beta, theta.hat_1se_halflambda = agg.re1_1se_halflambda$theta,
              beta.hat_min_lambda = agg.re1_min_lambda$beta, theta.hat_min_lambda = agg.re1_min_lambda$theta, 
              beta.hat_1se_lambda = agg.re1_1se_lambda$beta, theta.hat_1se_lambda = agg.re1_1se_lambda$theta,
              c1_lambda_const = c1_lambda_const))
}



## This performs Transfer Learning Lasso while estimating the Informative Auxiliary Data
## X is the matrix of covariates (nxp)
## y is the outcome vector (nx1)
## n.vec is a vector of number of observations in the target and aux datasets in that order. 
## RhatCount is either "n0/3" or a value between 1, ..., p to specify how many
## marginal correlations to consider when calculating information. 
## LambdaType is either "Constant" or "CV" to indicate if the lambda calculated 
## from the target dataset should be used and adjusted based on the aux dataset size ("Constant")
## or whether the optimal lambda should be calculated via cross validation ("CV") for
## each set of auxiliary information.
TransLasso.EstA0 <- function(X, y, n.vec, RhatCount, LambdaType, ...){
  
  # count of aux datasets
  M = length(n.vec)-1
  
  Rhat <- rep(0, M+1)
  p <- ncol(X)
  
  # make row indices of where target observations are
  ind.1 <- ind.set(n.vec, 1) 
  
  # calculate informativeness for each aux study
  for(k in 2: (M+1)){
    ind.k <- ind.set(n.vec, k) # row indices for kth aux sample. 
    
    # calculate difference in marginal correlations between k aux and target data. 
    Xty.k <- t(X[ind.k, ])%*%y[ind.k] / n.vec[k] - t(X[ind.1,])%*%y[ind.1]/ n.vec[1] 
    
    # change the number here to adjust how many marginal correlations are looked at!
    # take the top largest correlation differences based on number supplied in RhatCount
    margin.T <- sort(abs(Xty.k), decreasing=T)[1:RhatCount] 
    
    # estimated sparse index for kth aux sample.
    Rhat[k] <- sum(margin.T^2)  
  }
  
  Tset <- list()
  k0 = 0
  kk.list <- unique(rank(Rhat[-1])) # get ordering of smallest to largest Rhat for aux samples. 

  for(kk in 1:length(kk.list)){#use Rhat as the selection rule
    Tset[[k0+kk]] <- which(rank(Rhat[-1]) <= kk.list[kk])
  } # the sets of aux datasets to take for each ranking of datasets to include. 
  
  k0 = length(Tset)
  Tset <- unique(Tset)

  beta.T_min <- beta.T_min_lambda <- beta.T_min_halflambda <- list()
  beta.T_1se <- beta.T_1se_lambda <- beta.T_1se_halflambda <- list()
  
  # Lasso on Target Data only 
  init.re <- TL_Lasso(X=X, y=y, A0=NULL, n.vec=n.vec)
  
  beta.T_min[[1]] <- init.re$beta.kA_min 
  beta.T_1se[[1]] <- init.re$beta.kA_1se
  
  # if constant lambda specified, it is here. 
  c1_lambda_const <- init.re$lam.const
  c1_lambda_const_1se <- init.re$lam.const_1se
  
  og_lasso_coef_min <- init.re$beta.kA_min
  og_lasso_coef_1se <- init.re$beta.kA_1se
  
  beta.T_min_lambda <- beta.T_min_halflambda <- beta.T_min
  beta.T_1se_lambda <- beta.T_1se_halflambda <- beta.T_1se
  
  beta.pool.T_min <- beta.pool.T_min_lambda <- beta.pool.T_min_halflambda <- beta.T_min 
  beta.pool.T_1se <- beta.pool.T_1se_lambda <- beta.pool.T_1se_halflambda <- beta.T_1se 
  
  # go through TL Lasso for each informative set
  for(kk in 1:length(Tset)){ 
    T.k <- Tset[[kk]] 
    
    # which lambda type changes which section of function call it goes into
    if(LambdaType == "CV"){
      re.k <- TL_Lasso(X=X, y=y, A0=T.k, n.vec=n.vec, lam.const = NULL)
    }
    if(LambdaType == "Constant"){
      re.k <- TL_Lasso(X=X, y=y, A0=T.k, n.vec=n.vec, 
                          lam.const = c1_lambda_const, lam.const_1se = c1_lambda_const_1se)
    }
    
    # extract coefficients for each informative auxiliary set
    beta.T_min[[kk+1]] <- re.k$beta.kA_min 
    beta.pool.T_min[[kk+1]] <- re.k$w.kA_min
    
    beta.T_1se[[kk+1]] <- re.k$beta.kA_1se 
    beta.pool.T_1se[[kk+1]] <- re.k$w.kA_1se
    
    beta.T_min_halflambda[[kk+1]] <- re.k$beta.kA_min_halflambda 
    beta.pool.T_min_halflambda[[kk+1]] <- re.k$w.kA_min_halflambda
    
    beta.T_1se_lambda[[kk+1]] <- re.k$beta.kA_1se_lambda 
    beta.pool.T_1se_lambda[[kk+1]] <- re.k$w.kA_1se_lambda
  }
  
  beta.T_min <- beta.T_min[!duplicated((beta.T_min))] 
  beta.T_min <- as.matrix(as.data.frame(beta.T_min)) 
  
  beta.T_1se <- beta.T_1se[!duplicated((beta.T_1se))] 
  beta.T_1se <- as.matrix(as.data.frame(beta.T_1se)) 
  
  beta.T_min_halflambda <- beta.T_min_halflambda[!duplicated((beta.T_min_halflambda))] 
  beta.T_min_halflambda <- as.matrix(as.data.frame(beta.T_min_halflambda)) 
  beta.T_min_lambda <- beta.T_min_lambda[!duplicated((beta.T_min_lambda))] 
  beta.T_min_lambda <- as.matrix(as.data.frame(beta.T_min_lambda)) 
  
  beta.T_1se_halflambda <- beta.T_1se_halflambda[!duplicated((beta.T_1se_halflambda))] 
  beta.T_1se_halflambda <- as.matrix(as.data.frame(beta.T_1se_halflambda))
  beta.T_1se_lambda <- beta.T_1se_lambda[!duplicated((beta.T_1se_lambda))] 
  beta.T_1se_lambda <- as.matrix(as.data.frame(beta.T_1se_lambda)) 
  
  ## aggregate coefficients using squared error. 
  # No coef thresholding
  agg.re1_min <- coef.aggr(B= beta.T_min, X = X, y = y, N_vector = n.vec)
  agg.re1_1se <- coef.aggr(B= beta.T_1se, X = X, y = y, N_vector = n.vec)
  
  # half lambda coef thresholding 
  agg.re1_min_halflambda <- coef.aggr(B= beta.T_min_halflambda, X = X, y = y, N_vector = n.vec)
  agg.re1_1se_halflambda <- coef.aggr(B= beta.T_1se_halflambda, X = X, y = y, N_vector = n.vec)
  
  # lambda coef thresholding 
  agg.re1_min_lambda <- coef.aggr(B= beta.T_min_lambda, X = X, y = y, N_vector = n.vec)
  agg.re1_1se_lambda <- coef.aggr(B= beta.T_1se_lambda, X = X, y = y, N_vector = n.vec)
  
  
  # theta are the returned weights of each dataset. 
  # betas are the final, aggregated coefficients for potential covariates
  return(list(beta.hat_min = agg.re1_min$beta, theta.hat_min = agg.re1_min$theta, 
              beta.hat_1se = agg.re1_1se$beta, theta.hat_1se = agg.re1_1se$theta,
              beta.hat_min_halflambda = agg.re1_min_halflambda$beta, theta.hat_min_halflambda = agg.re1_min_halflambda$theta, 
              beta.hat_1se_halflambda = agg.re1_1se_halflambda$beta, theta.hat_1se_halflambda = agg.re1_1se_halflambda$theta,
              beta.hat_min_lambda = agg.re1_min_lambda$beta, theta.hat_min_lambda = agg.re1_min_lambda$theta, 
              beta.hat_1se_lambda = agg.re1_1se_lambda$beta, theta.hat_1se_lambda = agg.re1_1se_lambda$theta,
              c1_lambda_const = c1_lambda_const, 
              og_lasso_coef_min = og_lasso_coef_min, og_lasso_coef_1se = og_lasso_coef_1se))
}



################################################################################
############################ Complete Functions ################################


## This performs Transfer Learning Lasso allowing either the Auxiliary Dataset 
## Information to be known (Oracle or Oracle 1df) or estimated in the process 
## (Estimate A0). while estimating the Informative Auxiliary Data
## X_matrix is the matrix of covariates (nxp) concatenated between the Target 
## and all Auxiliary Datasets, in that order. n should therefore be n_0 + n_1 + ... + n_k.
## Please make sure that X_matrix is a matrix and not a dataframe to allow for 
## matrix multiplication to be carried out in the function. 
## Y_vector is the outcome vector to develop the TL Lasso model to estimate. It should 
## align with your X_matrix with Target and Auxiliary outcomes concatenated and be a nx1 vector.
## N_vector is a vector of number of observations in the Target and each Auxiliary
## datasets in the order they are concatenated. 
## AuxInformation is either "Estimate A0", "Oracle", or "Oracle 1df" to signify the 
## informative auxiliary datasets need estimated (EstA0) or they should be treated
## as known. If known, they can either be treated as equally informative and combined 
## into 1 dataset (Oracle), or they can be treated individually (Oracle 1df).
## RhatCount is either "n0/3" or a value between 1, ..., p to specify how many
## marginal correlations to consider when calculating information. This should be set to NULL
## if AuxInformation is "Oracle" or "Oracle 1df".
## LambdaType is either "Constant" or "CV" to indicate if the lambda calculated 
## from the target dataset should be used and adjusted based on the aux dataset size ("Constant")
## or whether the optimal lambda should be calculated via cross validation ("CV") for
## each set of auxiliary information. Default is "CV"
## seedstart is the seed to set in the calculation for reproducibility. If not specified,
## it is set to 123.

## This function will return a data.frame with TL coefficients in the columns. 
## The first column, "Variable" labels the intercept and column from your X matrix
## that it corresponds to. Columns 2:7 or 2:4 have the final coefficients with slight 
## variations in their calculation. "min" and "1se" correspond to the lambda that 
## either minimizes CV error or is at most 1se above it. "allcoef", "halfcoef",
## and "lambcoef" correspond to the coefficient thresholding used before 
## combining coefficients across the auxiliary and target dataset. If "Constant"
## LambdaType was specified, only the minimum lambda from the target data is used
## and so "1se" coefficients are not presented.
Epigen.TL.Lasso <- function(X_matrix, Y_vector, N_vector,
                           AuxInformation, RhatCount = NULL, LambdaType = "CV", 
                           seedstart = 123) {

  
  if(AuxInformation %!in% c("Estimate A0", "Oracle", "Oracle 1df")){
    stop("You must specify a method to capture auxiliary data informativeness.
         Options include 'Estimate A0', 'Oracle', or 'Oracle 1df'")
  }
  
  if(AuxInformation == "Estimate A0" & is.null(RhatCount)){
    stop("You must provide the number of Rhats to use when calculating auxiliary data informativeness.
         Options include any integer up to the column size of X_matrix OR 'n0/3'")
  }
  
  if(AuxInformation == "Estimate A0" & !is.null(RhatCount)){
    if(RhatCount > dim(X_matrix)[2]){
      stop("Rhat must be smaller than the number of columns in X. 
         Please specify any integer up to the column size of X_matrix OR 'n0/3'")
    }
    
    if(RhatCount == 0){
      stop("Rhat must be a positive number from 1, ..., column size of X_matrix.")
    }
  }
  
  
  if(LambdaType %!in% c("CV", "Constant")){
    stop("LambdaType must be either 'CV' or 'Constant' to specify which lambda parameter
         is used for the auxiliary data.")
  }
  
  if(is.null(X_matrix)){
    stop("Please supply the X matrix")
  } 
  if(is.null(Y_vector)){
    stop("Please supply the outcome")
  }
  if(is.null(N_vector)){
    stop("Please supply the N vector specifying the number of observations in target and auxiliary datasets")
  }
  
  if(sum(N_vector) != dim(X_matrix)[1]){
    stop("N_vector and X_matrix do not have the same number of observations")
  }
  
  if(sum(N_vector) != length(Y_vector)){
    stop("N_vector and Y do not have the same number of observations")
  }
  
  if(!is.null(RhatCount)){
    if(RhatCount == 'n0/3'){
        # set before calling the other functions. 
        RhatCount <- round(N_vector[1] / 3)
    }
  }
  
  
    # TransLasso for Estimating informative sex
    if(AuxInformation == "Estimate A0"){
      # for reproducibility
      set.seed(seedstart)

      Translasso_output <- TransLasso.EstA0(X = X_matrix, y = Y_vector,
                                    n.vec = N_vector, AuxInformation = AuxInformation,
                                          LambdaType = LambdaType, RhatCount = RhatCount)
    }
    # TransLasso for Oracle or Oracle 1df
    if(AuxInformation != "Estimate A0"){
      # for reproducibility
      set.seed(seedstart)
      
      Translasso_output <- TransLasso.Oracles(X = X_matrix, y = Y_vector, 
                                    n.vec = N_vector, AuxInformation = AuxInformation,
                                          LambdaType = LambdaType)
    }
    
    
    # calculate intercepts 
    intercept_beta_min <- mean(Y_vector - X_matrix %*% Translasso_output$beta.hat_min, na.rm=TRUE)
    intercept_beta_1se <- mean(Y_vector - X_matrix %*% Translasso_output$beta.hat_1se, na.rm=TRUE)

    intercept_beta_min_lambda <- mean(Y_vector - X_matrix %*% Translasso_output$beta.hat_min_lambda, na.rm=TRUE)
    intercept_beta_1se_lambda <- mean(Y_vector - X_matrix %*% Translasso_output$beta.hat_1se_lambda, na.rm=TRUE)
    
    intercept_beta_min_halflambda <- mean(Y_vector - X_matrix %*% Translasso_output$beta.hat_min_halflambda, na.rm=TRUE)
    intercept_beta_1se_halflambda <- mean(Y_vector - X_matrix %*% Translasso_output$beta.hat_1se_halflambda, na.rm=TRUE)
   
   
    # keep all coefficients. 
    TL_beta_coef <- data.frame(Variable = c("Intercept", colnames(X_matrix)),
                                    beta_min_allcoef = c(intercept_beta_min, Translasso_output$beta.hat_min),
                                    beta_1se_allcoef = c(intercept_beta_1se, Translasso_output$beta.hat_1se),
                                    beta_min_halfcoef = c(intercept_beta_min_halflambda, Translasso_output$beta.hat_min_halflambda),
                                    beta_1se_halfcoef = c(intercept_beta_1se_halflambda, Translasso_output$beta.hat_1se_halflambda),
                                    beta_min_lambcoef = c(intercept_beta_min_lambda, Translasso_output$beta.hat_min_lambda),
                                    beta_1se_lambcoef = c(intercept_beta_1se_lambda, Translasso_output$beta.hat_1se_lambda))
    
    # if Constant, only the min lambda is used. 
    if(LambdaType == "Constant"){
          TL_beta_coef <- data.frame(Variable = c("Intercept", colnames(X_matrix)),
                             beta_min_allcoef = c(intercept_beta_min, Translasso_output$beta.hat_min),
                             beta_min_halfcoef = c(intercept_beta_min_halflambda, Translasso_output$beta.hat_min_halflambda),
                             beta_min_lambcoef = c(intercept_beta_min_lambda, Translasso_output$beta.hat_min_lambda))
    } 
  
  # return Final TL coefficients as a dataframe
  return(TL_beta_coef)

}




## X is the matrix or dataframe of saliva DNA methylation beta values with 
## samples in rows and methylation sites in columns. To see the list of CpG 
## loci needed for computation, please call colnames to CS_Algorithms_GitHub or 
## C_Algorithms_GitHub. 
## method should be either "C+S" or "C" to indicate which set of algorithms 
## you are interested in calculating. If "C+S", you must also provide SalivaDNAmBiom. 
## SalivaDNAmBiom should be the matrix or dataframe of Saliva DNAm biomarkers 
## (ie DNAm biomarkers directly calculated with saliva methylation values) if
## using the C+S algorithms. Otherwise, this should be NULL. Default is NULL. 

## This function outputs a dataframe with predicted DNAm Biomarkers in each 
## column with rows in the same order as the originally supplied X. 

Saliva.2.Blood.DNAmBiomarkers <- function(X, method, SalivaDNAmBiom = NULL) {
  
  # make sure ppl provide the saliva DNAm biomarkers if its C+S
  if (method == "C+S" && is.null(SalivaDNAmBiom)) {
    stop("For method C+S, saliva DNAm Biomarker matrix must be provided.")
  }
  
  
  if (method == "C") {
    TLcoeffMatrix <- C_Algorithms_GitHub
    Varstart <- 2
  } else if (method == "C+S") {
    TLcoeffMatrix <- CS_Algorithms_GitHub
    Varstart <- 3
  } else {
    stop("Invalid method specified. Please specify either 'C' or 'C+S'")
  }
  
  
  cpgs_needed <- TLcoeffMatrix$Variable[Varstart:length(TLcoeffMatrix$Variable)]
  if(sum(cpgs_needed %!in% colnames(X)) > 0){
    cpgs_missing <<- cpgs_needed[cpgs_needed %!in% colnames(X)]
    stop("Not all CpG columns are present in X for this prediction. 
         Please see cpgs_missing for a list of missing but necessary CpGs.")
  }
  
   
  X_keep <- X[, cpgs_needed]
  if (any(is.na(X_keep))) {
    stop("Missing values are not allowed. Please impute missing values in the X matrix.")
  }
  
  ### Now can start the computations ###
  
  new_biomarker_list <- colnames(TLcoeffMatrix)[-1]
  
  n_biom <- dim(TLcoeffMatrix)[2]
  n_var <- dim(TLcoeffMatrix)[1]
  
  # make temp dataset with intercept, saliva, and columns in correct order
  if (method == "C"){
    newdata_temp <- data.frame(intercept = 1, X[, TLcoeffMatrix$Variable[Varstart:n_var]])
  }else{
    newdata_temp <- data.frame(intercept = 1, SalivaDNAmBiom = NA, 
                               X[, TLcoeffMatrix$Variable[Varstart:n_var]])
  }
  
  # Turn Coefficient matrix into a real matrix for multiplying
  TLcoeffMatrix2 <- matrix(unlist(TLcoeffMatrix[, 2:n_biom]), ncol = (n_biom-1), nrow = n_var)
  
  # initialize predictions dataframe 
  pred_df <- data.frame(rep(NA, dim(X)[1]))
  k <- 1
  
  for (i in new_biomarker_list) {
    
    if(method == "C+S"){
      # set Saliva value to the biomarker of interest using the name of the biomarker
      newdata_temp$SalivaDNAmBiom <- SalivaDNAmBiom[, i]
    }
    
    # get the column of the TL coef matrix the biomarker is in
    TL_coef_Col <- (which(colnames(TLcoeffMatrix) == i) - 1)
    
    # turn into matrix for multiplication 
    newdata_temp2 <- matrix(unlist(newdata_temp), ncol = dim(newdata_temp)[2], 
                            nrow = dim(newdata_temp)[1])
    
    # data has ppl in rows, variables in columns, 
    # TLcoef matrix has coefficients for each biomarker in a column
    cur_pred <- newdata_temp2 %*% TLcoeffMatrix2[, TL_coef_Col]
    
    # set prediction to the dataframe
    pred_df[, k] <- c(cur_pred)
    colnames(pred_df)[k] <- i 
    
    # make k go up 
    k <- k + 1
    
    # remove values so we don't repeat by accident
    rm(cur_pred); rm(newdata_temp2)
    
    if(method == "C+S"){
      # set Saliva value to NA so we don't repeat by accident
      newdata_temp$SalivaDNAmBiom <- NA
    }
    
  } # end of for loop in biomarker list 
  
  # return the predicted DNAm values. 
  if(method == "C+S"){
    colnames(pred_df) <- paste0(colnames(pred_df), "_CS_Pred")
  } else{
    colnames(pred_df) <- paste0(colnames(pred_df), "_C_Pred")
  }
  
  return(pred_df)
}



