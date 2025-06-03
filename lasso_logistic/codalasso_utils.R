
library(BenchmarkDenoise)
library(parallelly)
library(doParallel)
library(foreach)

repeat_codalasso <- function(y, X, lambdas=seq(0.05, 0.3, 0.05), train_prop=0.6, times=50, ncores=4){
  
  train_auc_mat <- matrix(0, nrow=length(lambdas), ncol=times)
  test_auc_mat <- matrix(0, nrow=length(lambdas), ncol=times)
  coefs_list <- list()
  
  # numCores <- availableCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  for(i in 1:length(lambdas)){
    current_lambda <- lambdas[i]
    print(current_lambda)
    coefs <- matrix(0, nrow=times, ncol=ncol(X))
    
    models <- foreach(j=1:times, .export=c("fit_codalasso")) %dopar%{
      
      fitted_model <- fit_codalasso(y=y, 
                               X=X,
                               lambdas=current_lambda,
                               train_prop=train_prop,
                               seed=j)
      fitted_model
      
    }
    
    for (j in 1:times){
      train_auc_mat[i, j] <- models[[j]]$train_auc
      test_auc_mat[i, j] <- models[[j]]$test_auc
      coefs[j, ] <- models[[j]]$coefs
    }
    coefs_list[[i]] <- coefs
  }
  
  stopCluster(cl)
  
  output <- list(lambdas=lambdas, 
                 train_auc=train_auc_mat,
                 test_auc=test_auc_mat,
                 coefs_list=coefs_list)
  
  return(output)
  
}

