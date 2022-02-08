\name{ProgPredLasso}
\alias{ProgPredLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identification of prognostic and predictive biomarkers
}
\description{
The computes the regularization path of the Prognostic Predictive Lasso described in the paper Zhu et al. (2022) given in the references.
}
\usage{
ProgPredLasso(X1, X2, Y=Y, cor_matrix=NULL, gamma=0.99, maxsteps=500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X1}{
Design matrix of patients characteristics with treatment 1
}
  \item{X2}{
Design matrix of patients characteristics with treatment 2
}
  \item{Y}{
Response variable
}
  \item{cor_matrix}{
Correlation matrix of biomarkers. If not specified, the function \code{\link{cvCovEst}} from package \code{\link{cvCovEst}} will be used to estimate this matrix.
}
  \item{gamma}{
Parameter \eqn{\gamma} defined in the paper Zhu et al. (2020) given in the references. Its default value is 0.99.}
  \item{maxsteps}{
Integer specifying the maximum number of steps for the generalized Lasso algorithm. Its default value is 500.}

}

\value{
Returns a list with the following components
\item{lambda }{different values of the parameter \eqn{\lambda} considered.}
\item{beta }{matrix of the estimations of \eqn{\beta} for all the \eqn{\lambda} considered.}
\item{beta.min }{estimation of \eqn{\beta} which minimize the MSE.}
\item{bic }{BIC for all the \eqn{\lambda} considered.}
\item{mse }{MSE for all the \eqn{\lambda} considered.}
}

\author{
Wencan Zhu, Celine Levy-Leduc, Nils Ternes
}

\examples{
X1 = t(sapply(c(1:25),FUN=function(x) rnorm(50)))
X2 = t(sapply(c(1:25),FUN=function(x) rnorm(50)))
Y=rnorm(50)
ProgPredLasso(X1=X1, X2=X2, Y=Y)

## The function is currently defined as
function(X1, X2, Y=Y, cor_matrix=NULL, gamma=0.99, maxsteps=500){
  p10=ncol(X1)
  p20=ncol(X2)
  p2 <- p10+p20
  
  mat_coef <- matrix(0,p2, p2)
  mat_coef[1:p10, 1:p10] <- diag(rep(1,p10))
  mat_coef[(p10+1):(p2), 1:p10] <- diag(rep(-1,p10))
  mat_coef[(p10+1):(p2), (p10+1):(p2)] <- diag(rep(1,p10))
  
  n1=nrow(X1)
  n2=nrow(X2)
  n=n1+n2
  
  if((n1+n2)!=length(Y))
    stop("the sample size should be consistent")
  
  TRT1  <- c(rep(1,n1), rep(0, n2))
  TRT2  <- c(rep(0,n1), rep(1, n2))
  
  X_full <- cbind(rbind(X1, X2)*TRT1, rbind(X1, X2)*TRT2)
  X_classic <- cbind(rbind(X1, X2), rbind(X1, X2)*TRT2)
  
  if(is.null(cor_matrix)){
    X_all <- rbind(X1, X2)
    cv_cov_est_out <- cvCovEst(
      dat = X_all,
      estimators = c(
        linearShrinkLWEst, denseLinearShrinkEst,
        thresholdingEst, poetEst, sampleCovEst
      ),
      estimator_params = list(
        thresholdingEst = list(gamma = c(0.2, 0.4)),
        poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L))
      ),
      cv_loss = cvMatrixFrobeniusLoss,
      cv_scheme = "v_fold",
      v_folds = 5
    )
    cor_matrix <- cov2cor(cv_cov_est_out$estimate)
  }
  
  cor_matrix_full = matrix(0, (p2), (p2))
  cor_matrix_full[1:p10, 1:p10] = cor_matrix_full[((p10+1):(p2)), ((p10+1):(p2))] = cor_matrix
  
  cor_matrix_full <- round(cor_matrix_full, 6)
  file <- try(SVD_sigma <- svd(cor_matrix_full))
  if (class(file) == "try-error") {
    cat("Caught an error during SVD.\n")
    Eigen_Sigma <- eigen(cor_matrix_full)
    V_sigma <- Eigen_Sigma$vectors
    lam <- Eigen_Sigma$values
    square_root_sigma <- V_sigma%*%diag(sqrt(lam))%*%solve(V_sigma)
    inv_diag <- ifelse(lam<0.000001, 0, 1/sqrt(lam))
    inv_square_root_Sigma <- V_sigma%*%diag(inv_diag)%*%solve(V_sigma)
  } else {
    U_sigma <- SVD_sigma$u
    D_sigma <- SVD_sigma$d
    square_root_sigma <- U_sigma%*%diag(sqrt(D_sigma))%*%t(U_sigma)
    inv_diag <- ifelse(D_sigma<0.000001, 0, 1/sqrt(D_sigma))
    inv_square_root_Sigma <- U_sigma%*%diag(inv_diag)%*%t(U_sigma)
    inv_square_root_Sigma <- U_sigma%*%diag(1/sqrt(D_sigma))%*%t(U_sigma)
  }

  #tranformation matrix
  inv_square_root_Sigma_trt = diag(1, ncol=(p2+2), nrow=(p2+2))
  inv_square_root_Sigma_trt[3:(p2+2), 3:(p2+2)] <- inv_square_root_Sigma
  
  mat_trt <- matrix(0,nrow=p2, ncol=(p2+2) )
  mat_trt[, 3:(p2+2)] <- mat_coef
  mat_trt[, c(1,2)] <- matrix(0, p2, 2)
  
  X_new <- cbind(c(rep(1, n1),rep(0, n2)), c(rep(0, n1),rep(1, n2)), X_full)
  X_new2 <- cbind(c(rep(1, n1),rep(0, n2)), c(rep(0, n1),rep(1, n2)), X_classic)
  trans_mat <- mat_trt\%*\%inv_square_root_Sigma_trt
  X_tilde0 <- X_new\%*\%inv_square_root_Sigma_trt 
  out0 <-  genlasso(Y, X_tilde0, trans_mat, maxsteps=maxsteps)

 if(p10<=50){ 
    top_grill <- seq(1, p10, 2)
  } else if(p10<210){ 
    top_grill <- c(1:50, seq(52,p10, 2))
  } else if (p10<500){
    top_grill <- c(1:50, seq(52,100, 2), seq(105,200, 5), seq(210, p10, 10))
  } else {
    top_grill <- c(1:50, seq(52,100, 2), seq(105,200, 5), seq(210, 500, 10))
  }

  opt_top <- opt_final_top <- c()
  beta_final <- matrix(NA, length(out0$lambda), (p2+2))
  mse_final = bic_final = c()
  for(i in 1:length(out0$lambda)){
    #### Filter on beta_tilde
    gamma_tilde <- out0$beta[,i][-c(1,2)]
    gamma_tilde_prog <- gamma_tilde[1:p10]
    gamma_tilde_pred <- gamma_tilde[(p10+1):(p2)]
    
    gamma_tilde_opt <-Correction2Vect(X=X_tilde0, Y=Y, te=out0$beta[c(1,2),i], 
    vector_prog =gamma_tilde_prog, vector_pred= gamma_tilde_pred ,delta = 0.99999, 
    top_grill. = top_grill)
    
    #### Filter on beta
    beta_final0 <- mat_coef\%*\%inv_square_root_Sigma\%*\%gamma_tilde_opt
    beta_prog <- beta_final0[1:p10]
    beta_pred <- beta_final0[(p10+1):(p2)]
    
    beta_opt_final <- Correction2Vect(X=X_new2, Y=Y, te=out0$beta[c(1,2),i], 
    vector_prog =beta_prog, vector_pred= beta_pred ,delta = 0.99999, 
    top_grill. = top_grill, toZero = TRUE)
    
    beta_final[i, ] <- round(c(out0$beta[c(1,2),i], beta_opt_final),6)
    
    X_temp <- X_classic[, which(beta_opt_final!=0)]
    if(length(which(beta_opt_final!=0))>=(length(Y)-2)){
      X_pred <- as.matrix(cbind(c(rep(1, n1),rep(0, n2)), c(rep(0, n1),rep(1, n2)), X_temp))
      mod_ridge <- cv.glmnet(x=as.matrix(X_pred), y=Y, alpha=0, intercept=FALSE)
      mse_final[i] <- min(mod_ridge$cvm)
      bic_final[i] <- n*log(mse_final[i])+log(n)*length(which(beta_opt_final!=0))
    } else {
      mydata <- data.frame(Y=Y, cbind(c(rep(1, n1),rep(0, n2)), c(rep(0, n1),rep(1, n2)), X_temp))
      formula <- paste0("Y~-1 +",paste0(colnames(mydata[, -which(colnames(mydata)=="Y")]), 
      collapse=" + "))
      myform <- as.formula(formula)
      mod_lm <- lm(myform, data=mydata)
      mse_final[i] <- mean(mod_lm$residuals^2)
      bic_final[i] <- n*log(mse_final[i])+log(n)*length(which(beta_opt_final!=0))
    }
  }
  
  beta_min <- beta_final[which.min(bic_final), ]
  return(list(lambda=out0$lambda, beta=beta_final, beta.min=beta_min, bic=bic_final, mse=mse_final))
}
}