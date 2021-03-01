library(mvtnorm)
library(parallel)
library(Matching)
library(regpro)
library(FNN)
library(ggplot2)
library(splines2)
library(base)
library(dplyr)
library(magrittr)
library(gam)
library(expm)
#Auxiliary Functions
#Smoothing splines 
fit.SS <- function(X, Y, pred_X){
  data_train <- data.frame(Y = Y, X = X)
  form<-as.formula(
    paste0("Y~",paste0("s(X.",seq( 1, ncol(X), 1),")",collapse="+")))
  model_ET_X <- gam(data=data_train, formula=form)
  predict(object = model_ET_X,type="response",newdata = data.frame(X = pred_X)) %>% as.vector()
}
##KNN regression and CV
make_knn_pred = function(k, training, testing, y_ta, y_tt) {
  pred = FNN::knn.reg(train = training, 
                      test = testing, 
                      y = y_ta, k = k)$pred
  return(mean((pred - y_tt)^2))
}
CV_knn <- function(j,k_pool,x,y){
  error_k <- vector()
  for (k in k_pool){
  n <- length(y)
  nset <- seq(1,n,1)
  nj <- floor(n/j)
  subset <- list()
  for(i in 1:(j-1)) {
    subset[[i]] <- sample(nset,nj)
    nset <- nset[-which(nset %in% subset[[i]])]
  }
  subset[[j]] <- nset
  nset <- seq(1,n,1)
  error <- vector()
  for (i in 1:j){
    test <- subset[[i]]
    train <- nset[-which(nset %in% subset[[i]])]
    train_x <- x[train,]; test_x = x[test,]
    train_y <- y[train]; test_y = y[test]
    error <- c(error,make_knn_pred(k,train_x,test_x,train_y,test_y))
  }
  error_k <- c(error_k, mean(error))
  }
  return(k_pool[which(error_k == min(error_k))][1])
}
##Basis construction
Bspline_construct <- function(X ,df, degree, range){
  knots <- seq(range[1], range[2], (range[2] - range[1])/(df - degree))
  knots <- knots[-1]; knots <- knots[1:(df - degree-1)] 
  return(bSpline(X, knots = knots, degree = degree, intercept = TRUE, Boundary.knots = range))
}
Powerseries_construct_T = function(X ,df, range){
  mid_point <- (range[1] + range[2])/2
  range <- (range[2] - range[1])/2
  M <- matrix(0,ncol = df, nrow = length(X))
  for (i in 1:length(X)){
    for (j in 1:df){
      M[i,j] <- (X[i] - mid_point)^j/(2*(range)^(j + 1)/(j + 1))
    }
  }
  return(M)
}
Powerseries_construct_V = function(X ,df, range){
  mid_point <- (range[1] + range[2])/2
  range <- (range[2] - range[1])/2
  M <- matrix(0,ncol = df, nrow = length(X))
  for (i in 1:length(X)){
    for (j in 1:df){
      M[i,j] <- (X[i] - mid_point)^(j - 1)/(2*(range)^(j)/(j))
    }
  }
  return(M)
}
rowwise_kron <- function(M1,M2){
  M <- matrix(0, nrow <- nrow(M1), ncol <- ncol(M1)*ncol(M2))
  for (i in 1:nrow(M1)){
    M[i,] <- kronecker(M1[i,],M2[i,])
  }
  return(M)
}
Basis_M_Build <- function(V, Tt, Vbasis_name, Tbasis_name, df_V, df_T, degree, Vrange, Trange){
  n <- length(Tt); dim_V <- ncol(V)
  basis_V <-  case_when(
    Vbasis_name == 'Bspline' ~ {
      if (Vbasis_name == 'Bspline' ){
      apply(V, MARGIN = 2, FUN = function(v){Bspline_construct(v, df_V, degree, Vrange)})
      } else {0}
    },
    Vbasis_name == 'PowerSeries' ~ {
      apply(V, MARGIN = 2, FUN = function(v){Powerseries_construct_V(v, df_V, Vrange)})
    }
  ) %>% 
    as.data.frame %>% as.list %>%
    lapply(function(v){v %>% as.vector %>% matrix(nrow = n)})
  
  basis_T <- case_when(
    Tbasis_name == 'Bspline' ~ {
      # apply(Tt, MARGIN = 2, FUN = function(v){Bspline_construct(v, df_T, degree, Trange)})
      0
    },
    Tbasis_name == 'PowerSeries' ~ {
      apply(Tt, MARGIN = 2, FUN = function(v){Powerseries_construct_T(v, df_T, Trange)})
    }
  ) %>% 
    as.data.frame %>% as.list %>%
    lapply(function(v){v %>% as.vector %>% matrix(nrow = n)})
  
  basis_M <- matrix(1, nrow = n, ncol = 1)
  for (i in 1:dim_V){
    basis_M <- rowwise_kron(basis_M, basis_V[[i]])
  }
  basis_M <- rowwise_kron(basis_M, basis_T[[1]])
  return(basis_M)
}


#Main Functions
GCCTE_fit = function(X, idxV = c(1), Tt, Y, Vbasis_name ='Bspline', Tbasis_name = 'PowerSeries',
                     df_V_pool, df_T_pool, degree, q_type = 'optimal', mu0, ET_X, Vrange, Trange, k_pool, J = 7, 
                     nusi_est = 'SS', CI = FALSE){
#variable description
##X: n by dx matrix; idxV: coordinates in X that are V; Tt: treatment vector; Y: outcome vector
##Vbasis: basis type for V; Tbasis:basis type for T
##df_V_pool, df_T_pool: degrees of V that one has interest to test with GCV 
##degree: degree of the Bspline
##Kn_max: max number of Kn for GCV
##q_type: type of q function
##mu0, ET_X: estimated nusicance vectors
##Vrange, Trange: lists for range of V and T
##k_pool: candidates of k for knn regression
##J: number of rounds for knn's CV
##nusi_est: estimator for nusicance, 1) 'SS' = smoothness spline with default settings in  GAM package; 2) 'KNN' = k-nearest neighborhood regression
##CI: TRUE if want to get confidence interval in est_tau
  
#solve function to fit with certian degree of V and T
  solve_GCCTE = function(df_V, df_T, CI = FALSE){
    Kn = df_V * df_T
    #build basis
    basis_M <- Basis_M_Build(V, Tt, Vbasis_name, Tbasis_name, df_V, df_T, degree, Vrange, Trange)
    #build Q matrix of q function
    Q <- case_when( 
      q_type == 'optimal' ~ {
        q_train <- apply(basis_M, MARGIN = 2, function(v){(ET_X - Tt)*v}) %>%
          as.data.frame %>% as.list;
        if (nusi_est == 'KNN'){
        kopt_Q <- lapply(q_train, function(v){CV_knn(J, k_pool, X, as.matrix(v))});
        #combine the selected k and sample's q for each basis dimenison, to one list
        q_train_with_k <- Map(list, kopt_Q, q_train);
        Q_num <- q_train_with_k %>% lapply(FUN = function(l){FNN::knn.reg(train = X, 
                                                                          test = X, 
                                                                          y = l[[2]], k = l[[1]])$pred}) %>% unlist;
        }
        if (nusi_est == 'SS'){
          Q_num <- q_train %>% lapply(FUN = function(l){fit.SS(X, l, X)}) %>% unlist;
        }
        Q_num
      } 
    ) %>% matrix(nrow = n)
    #Get Smoothing matrix S for GCV and Phi for tau estimation
    QKn <- Q[, 1:Kn]
    MpartI <- apply(QKn, 2, function(v){v*(Tt - ET_X)}) %>% t %*% basis_M[, 1:Kn] %>% solve 
    MpartII <- QKn %>% t %*% diag(as.numeric(Tt - ET_X))
    MpartIII <- basis_M[, 1:Kn] 
    phi <- MpartI %*% MpartII %*% (Y - mu0)
    S <- MpartIII %*% MpartI %*% MpartII
    MSE <- mean((Y - mu0 - S %*% (Y - mu0))^2)
    #pre-A for CI
    if (CI){
      tau_samples <- S %*% (Y - mu0) %>% as.vector
      G = matrix(0,Kn,Kn)
      B = matrix(0,Kn,Kn)
      for (obj in 1:n){
        q_hat = Q[obj,] %>% as.matrix
        G = G + (Tt[obj] - ET_X[obj]) * q_hat %*% basis_M[obj,]
        B = B + (Y[obj] - tau_samples[obj] - mu0[obj])^2 * (Tt[obj] - ET_X[obj])^2 * q_hat %*% t(q_hat)
      }
      G = G/n
      B = B/n
      pre_A = 1/sqrt(Kn) * solve(G) %*% sqrtm(B)
      return(list(phi = phi, S = S, MSE = MSE, pre_A = pre_A))
    }else{
    return(list(phi = phi, S = S, MSE = MSE))
    }
  }
#pre-setting 
  V <- X[,idxV] %>% 
    as.matrix()
  dim_V <- length(idxV)
  n <- length(Y)
#GCV df selection
  error_GCV = matrix(0, nrow = length(df_V_pool), ncol = length(df_T_pool))
  for (df_V in df_V_pool){
    for (df_T in df_T_pool){
      solve_GCV <- solve_GCCTE(df_V, df_T, FALSE)
      S <- solve_GCV$S; phi = solve_GCV$phi 
      error_GCV[which(df_V_pool == df_V), which(df_T_pool == df_T)] <- solve_GCV$MSE/((1 - sum(diag(S))/n)^2)
    }
  }
  df_V_opt <- df_V_pool[which(error_GCV == min(error_GCV), arr.ind = TRUE)[1,1]]
  df_T_opt <- df_T_pool[which(error_GCV == min(error_GCV), arr.ind = TRUE)[1,2]]
#fitting
  final_solve <- solve_GCCTE(df_V_opt, df_T_opt, TRUE)
  final_model <- list(para = final_solve$phi, df_V = df_V_opt, df_T = df_T_opt, Vbasis_name = Vbasis_name,
                      Tbasis_name = Tbasis_name, degree = degree, Vrange = Vrange, Trange = Trange, pre_A = final_solve$pre_A)
  return(final_model)
}

est_tau <- function(V_pred, Tt_pred, model, CI = FALSE, rate = 0.9){
  basis_M_pred <- Basis_M_Build(V_pred, Tt_pred, model$Vbasis_name, model$Tbasis_name, model$df_V, 
                                model$df_T, model$degree, model$Vrange, model$Trange)
  tau = basis_M_pred %*% model$para
  z_value = qnorm(1 - (1 - rate)/2)
  if (CI){
    Kn <- model$df_V * model$df_T
    upper <- rep(0,length(Tt_pred)); lower <- rep(0,length(Tt_pred))
    for (obj in 1:length(Tt_pred)){
      A <- basis_M_pred[obj,] %*% model$pre_A
      est_gamma <- A %*% t(A)
      upper[obj] <- tau[obj] + z_value * sqrt(est_gamma * (Kn) / n)
      lower[obj] <- tau[obj] - z_value * sqrt(est_gamma * (Kn) / n)
    }
    return(list(tau = tau, upper = upper, lower = lower))
  }else{
  return(list(tau = tau))
  }
}

