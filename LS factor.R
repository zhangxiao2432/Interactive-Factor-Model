#====================================================================================
# LEAST SQUARES ESTIMATION OF LINEAR PANEL DATA MODELS WITH INTERACTIVE FIXED EFFECTS
#====================================================================================
#   
# AUTHOR OF THIS CODE: 
#=====================
#
#    Martin Weidner, University College London, and CeMMaP
#    Email: m.weidner (at) ucl (dot) ac (dot) uk
#    (last update: December 2013)
#   
#    DISCLAIMER:
#   ============
#   
#    This code is offered with no guarantees. Not all features of this code
#    were properly tested. Please let me know if you find any bugs or
#    encounter any problems while using this code. All feedback is
#    appreciated.
#   
#    REFERENCES:
#   ============
#   
#    For a description of the model and the least squares estimator see e.g.
#    Bai (2009, "Panel data models with interactive fixed effects"), or
#    Moon and Weidner (two working papers:
#                         % "Dynamic Linear Panel Regression Models with Interactive Fixed Effects"
#                       % "Linear Regression for Panel with Unknown Number of Factors as Interactive Fixed Effects").
#   
#    THREE DIFFERENT COMPUTATION METHODS ARE IMPLEMENTED:
#   =====================================================
#       
#     METHOD 1: (recommended default method)
#     ----------
#     iterate the following two steps until convergence:
#     Step 1: for given beta compute update for lambda and f as
#               principal components of Y-beta*X
#     Step 2: for given lambda and f update beta by runing a pooled OLS regression of
#               M_lambda*Y*M_f (or equivalently of just Y itself) on M_lambda*X*M_f.
#     The procedure is repeated multiple times with different starting values.
#     
#     
#      METHOD 2: 
#     ----------
#     The profile objective function (after profiling out lambda and f) is
#     optimized over beta using "fminunc"
#     The procedure is repeated multiple times with different starting values.
#     
#     
#     METHOD 3: (described in Bai, 2009)
#     ----------
#     iterate the following two steps until convergence:
#     Step 1: for given beta compute update for lambda and f as
#               principal components of Y-beta*X (same as in method 1)
#     Step 2: for given lambda and f run a pooled OLS regression of
#               Y-lambda*f' on X to update beta.
#     The procedure is repeated multiple times with different starting values.
#     
#     
#     COMMENTS:
#     ==========
#     Another method would be to use Step 1 as in Method 1&3, but to 
#     replace step 2 with a regression of Y on either M_lambda*X or X*M_f, i.e.
#     to only project out either lambda or f in the step 2 regression.
#     Bai (2009) mentions this method and refers to Ahn, Lee, and Schmidt (2001),
#     Kiefer (1980) and Sargan (1964) for this. We have not tested this
#     alternative method, but we suspect that Method 1 performs better in
#     terms of speed of convergence.
#     
#     This alternative method and the method proposed by Bai (2009) --- i.e.
#    "method 3" here --- have the property of reducing the LS objective function in
#     each step. This is not true for Method 1 and may be seen as a
#     disadvantage of Method 1. However, we found this to be a nice feature,
#     because we use this property of Method 1 as a stopping rule:
#     if the LS objective function does not improve, then we know we are
#    "far away" from a proper minimum, so we stop the iteration
#     and begin the iteration with another randomly chosen starting value.
#     Note that multiple runs with different starting values are required
#     anyways for all methods (because the LS objective function may have
#     multiple local minima).
#    
#     We recommend method 1, because each iteration step is fast (quicker
#     than method 2, which needs to compute a gradient (and Hessian?) in each
#     step, involving multiple evaluations of the objective function) and its
#     rate of convergence in our tests was very good (faster than method 3).
#     However, we have not much explored the relative sensativity of the
#     different methods towards the choice of starting value. Note that by
#     choosing the quickest method (method 1) one can try out more different
#     starting values of the procedure in the same amount of time.
#     Nevertheless, it may well be that method 2 or 3 or the alternative
#     method described above perform better in certain situations.

LS.factor <- function(Y,X,R,report=TRUE,precision_beta=10^-8,method="m1",start=matrix(0,dim(X)[3],1),repMIN=30,repMAX=10*repMIN,M1=1,M2=0){
  
  #============================================================================
  # REQUIRED INPUT PARAMETERS: 
  #     Y = NxT matrix of outcomes
  #     X = KxNxT multi-matrix of regressors
  #     R = positive integer 
  #         ... number of interactive fixed effects in the estimation
  #
  #     Comment: we assume a balanced panel, i.e. all elements of Y and X are known.
  #
  # OPTIONAL INPUT PARAMETERS: 
  # report = 'silent' ... the program is running silently
  #        = 'report' ... the program reports what it is doing
  # precision_beta = defines stopping criteria for numerical optimization, namely optimization is stopped when difference in beta
  #                  relative to previous opimtization step is smaller than "precision_beta" (uniformly over all K components of beta)
  #                  NOTE: the actual precision in beta will typically be lower than precision_beta, depending on the convergence rate of the procedure.
  #     
  #     method = 'm1' or 'm2' or 'm3' ... which optimization method is used 
  #                              (described above)
  #     start = Kx1 vector, first starting value for numerical optimization
  #     repMIN = positive integer
  #              ... minimal number of runs of optimization with different starting point
  #     repMAX = positive integer
  #              ... maximal number of runs of optimization (in case numerical optimization doesn't terminate properly, we do multiple runs even for repMIN=1)
  #     M1 = positive integer 
  #          ... bandwidth for bias correction for dynamic bias (bcorr1),
  #          M1=number of lags of correlation between regressors and errors that is corrected for in dynamic bias correction
  #     M2 = non-negative integer 
  #          ... bandwidth for bias correction for time-serial correlation (bcorr3),
  #          M2=0 only corrects for time-series heteroscedasticity, while M2>0 corrects for time-correlation in erros up to lag M2 
  #
  #     OUTPUT PARAMETERS:
  #     beta = parameter estimate 
  #     exitflag = 1 if iteration algorithm properly converged at optimal beta
  #              = -1 if iteration algorithm did not properly converge at optimal beta
  #     lambda= estimate for factor loading
  #     f = estimate for factors
  #     Vbeta1 = estimated variance-covariance matrix of beta,
  #              assuming homoscedasticity of errors in both dimensions
  #     Vbeta2 = estimated variance-covariance matrix of beta,
  #              assuming heteroscedasticity of errors in both dimensions
  #     Vbeta3 = estimated variance-covariance matrix of beta,
  #              allowing for time-serial correlation up to lag M2
  #              (i.e. if M2==0, then Vbeta2==Vbeta3)
  #     bcorr1,2,3 = estimate for the three different bias components
  #              (needs to be subtracted from beta to correct for the bias)
  #              bcorr1 = bias due to pre-determined regressors
  #              bcorr2 = bias due to cross-sectional heteroscedasticity of errors
  #              bcorr3 = bias due to time-serial heteroscedasticity and time-serial correlation of errors
  #============================================================================================================                                                               
  
  # Y <- as.matrix(Y[,-c(1:2)])
  # X <- Ylag.3d
  # R <- 1
  # report=TRUE
  # precision_beta=1e-4
  # method="m1"
  # start=matrix(0,dim(X)[3],1)
  # repMIN=30
  # repMAX=10*repMIN
  # M1=1
  # M2=0
  
  #COMMENT: We assume that all provided input parameters have values and dimensions as descibed above. The program could 
  #be impoved by checking that this is indeed the case.
                                                             
  K <- dim(X)[3]     #number of regressors
  N <- dim(X)[1]     #cross-sectional dimension
  TT <- dim(X)[2]     #time-serial dimension  
  
  #Input parameters that are not provided are given default parameters as follows:
  
  #if N<T we permute N and T in order to simplify computation of eigenvectors and eigenvalues (we only solve the 
  #eigenvalue problems for TxT matrices and we make sure that T<=N)
  
  trans <- 0
  if(N<TT){
    
    trans <- 1              #dummy variable to remember that we exchanged N and T dims.
    NN <- N
    N <- TT
    TT <- NN
    Y <- t(Y)
    X <- aperm(X,c(3,2,1))     
    
  }
  
  #NUMERICAL OPTIMIZATION TO OBTAIN beta:
    
  beta <- 1e5*matrix(1,length(start),1)          #best beta found so far
  obj0 <- 1e4                                    #objective function at optimal beta found so far
  count <- 0                                     #how many minimization runs properly converged so far
  exitflag <- -1                                 #no proper solution found, yet 
  
  for(i in 1:repMAX){
    if(count < repMIN){
      
      #CHOOSE STARTING VALUE FOR OPTIMIZATION
      if(i==1){
        st <- start              #first starting value is given by user (or =0 by default)
      } else {
        st <- start + 10*matrix(rnorm(length(start)),length(start),1) #choose random starting values up from second run
        #COMMENT: this is a very simple way of choosing different starting values. One might want 
        #to modify this deopending on the problem one consideres.
      }
      
      #REPORT TO USER:
      if(report==TRUE){
        
        cat("Program LS_factor now starting optimization run number: ", toString(i), "/", toString(repMAX),"\n")
        cat("number of optimization runs that converged so far: ", toString(count), "/", toString(repMIN),"\n")
        cat("starting value for current run = ", toString(st),"\n")
        
      }
      
      #RUN ACTUAL OPTIMIZATION OVER beta:
      if(method=="m1"){
        
        temp <- minimize_obj_method1(Y,X,R,st,precision_beta)
        para <- temp$beta
        obj <- temp$obj
        ef <- temp$ef
        iter <- temp$iter
        #cat("Iteration: ", iter, "\n")
        
      } else if(method=="m2") {
        
        temp <- minimize_obj_method2(Y,X,R,st,precision_beta)
        para <- temp$beta
        obj <- temp$obj
        
      } else {
        
        temp <- minimize_obj_method3(Y,X,R,st,precision_beta)
        para <- temp$beta
        obj <- temp$obj
        ef <- temp$ef
        iter <- temp$iter
        #cat("Iteration: ", iter, "\n")
        
      }
      
      #REPORT TO USER:
      if(report==TRUE){
        
        if(ef > 0){
          cat("Method", toString(method), "converged at beta = ", toString(t(para)),"\n")
        } else {
          cat("Method", toString(method), "did NOT converge. Stopped at beta = ", toString(t(para)),"\n")
        }
        
        if(obj < obj0){
          cat("Final Objective = ", toString(obj), " ==> NEW BEST MINIMUM FOUND", "\n")
        } else {
          cat("Final Objective = ", toString(obj), " > Best Objective so Far = ", toString(obj0),"\n")
        }
        
      }
      
      #UPDATE ESTIMATOR, IN CASE BETTER SOLUTION FOUND:
      if(obj < obj0){
        obj0 <- obj
        beta <- para       #new "global minimum" found
        
        if(ef > 0){
          exitflag <- 1       #optimal beta corresponds to point where iteration algorithm properly converged 
        } else {
          exitflag <- -1
        }
      }
      
      #UPDATE COUNTER OF "good" SOLUTIONS:
      if(ef > 0){                #%if method properly converged, then count how many "good" solutions found
        count <- count+1
      }
      
    }
  }          #%end of calculation of beta-estimator
  
  #CALCULATE lambda AND f FOR THE OPTIMAL beta:
  res1 <- Y
  
  for(k in 1:K){
    res1 <- res1 - beta[k]*drop(X[k,,])
  }
  
  eig <- eigen(t(res1)%*%res1)
  V <- eig$vectors
  D <- eig$values
  
  for(r in 1:R){
    
    V[,r] <- V[,r]/norm(as.matrix(V[,r]),type = "F")
    
    if(mean(V[,r]) < 0){
      V[,r] <- -V[,r]
    }
  }
  
  f <- V[,1:R]
  lambda <- res1%*%f
  res <- res1 - lambda%*%t(f)               #estimate for the residuals
  
  if(trans==1){                           #need to undo the interchange of N and T now
    
    save <- lambda
    lambda <- V
    V <- save    
    res <- t(res)
    NN <- N
    N <- TT
    TT <- NN
    Y <- t(Y)
    X <- aperm(X,c(3,2,1))     
    
  }
  
  #CALCULATE VARIANCE-COVARIANCE MATRIX OF beta:
  Pf <- f%*%solve(t(f)%*%f)%*%t(f)
  Plambda <- lambda%*%solve(t(lambda)%*%lambda)%*%t(lambda)
  Mf <- diag(TT) - Pf
  Mlambda <- diag(N) - Plambda
  W <- matrix(0,K,K)
  Omega <- Omega2 <- matrix(0,K,K)
  
  for (k1 in 1:K) {
    for (k2 in 1:K){
      
      Xk1 <- Mlambda%*%drop(X[,,k1])%*%Mf
      Xk2 <- Mlambda%*%drop(X[,,k2])%*%Mf
      W[k1,k2] <- 1/N/TT*trace1(Mf%*%t(drop(X[,,k1]))%*%Mlambda%*%drop(X[,,k2])) #Hessian
      Omega[k1,k2] <- 1/N/TT*t(c(Xk1)*c(Xk2))%*%(c(res)^2)  #Variance of Score
      Omega2[k1,k2] <- 1/N/TT*trace1(trunc1(t(res*Xk1)%*%(res*Xk1),M2+1,M2+1))
      
    }
  }
  
  sigma2 <- trace1(t(res)%*%res)/N/TT
  Vbeta1 <- solve(W)%*%sigma2/N/TT
  Vbeta2 <- solve(W)%*%Omega%*%solve(W)/N/TT
  Vbeta3 <- solve(W)%*%Omega2%*%solve(W)/N/TT
  B1 <- B2 <- B3 <- matrix(NA,K,1)
  
  for(k in 1:K){
    
    XX <- drop(X[,,k])
    B1[k] <- 1/sqrt(N*TT)%*%trace1(Pf%*%trunc1(t(res)%*%XX,M1,M1))
    B2[k] <- 1/sqrt(N*TT)%*%trace1(t(XX)%*%Mlambda%*%trunc1(res%*%t(res),1,1)%*%lambda%*%solve(t(lambda)%*%lambda)%*%solve(t(f)%*%f)%*%t(f))
    B3[k] <- 1/sqrt(N*TT)%*%trace1(trunc1(t(res)%*%res,M2+1,M2+1)%*%Mf%*%t(XX)%*%lambda%*%solve(t(lambda)%*%lambda)%*%solve(t(f)%*%f)%*%t(f))
    
  }
  
  den <- matrix(0,K,1)
  for(k1 in 1:K){
    for (k2 in 1:K) {
      
     den[k1,k2] <- 1/N/TT*trace1(Mf%*%t(drop(X[,,k1]))%*%Mlambda%*%drop(X[,,k2]))
      
    }
  }
  
  bcorr1 <- -solve(den)%*%t(B1)/sqrt(N*TT)
  bcorr2 <- -solve(den)%*%t(B2)/sqrt(N*TT)
  bcorr3 <- -solve(den)%*%t(B3)/sqrt(N*TT)
  
  return(list(beta=beta,exitflag=exitflag,lambda=lambda,f=f,Vbeta1=Vbeta1,Vbeta2=Vbeta2,Vbeta3=Vbeta3,bcorr1=bcorr1,bcorr2=bcorr2,bcorr3=bcorr3))
}


minimize_obj_method1 <- function(Y,X,R,st,precision_beta){
  
  #INPUT: Y  = NxT
  #       X  = KxNxT
  #       st = Kx1 ... starting value for optimization over regression paramter beta
  #       precision_beta = defines stopping criteria, namely optimization is stopped when difference in beta after one
  #                        optimization step is smaller than precision_beta (uniformly over all K components of beta)
  #OUTPUT: beta = Kx1 ... optimal beta that was found
  #        obj  = value of LS-objective function at optimal beta
  #        ef (exitflag) = 1  if procedure properly terminated (according to "precision_beta" criteria)
  #                      = -1 if procedure failed to converge (namely if
  #objective function did not improve during last step --- in principle the procedure could still converge afterwards, but if the
  #objective function does not improve, then it seems more promising to stop the optimization and restart with a different starting value)
  
  N <- dim(Y)[1]
  TT <- dim(Y)[2]
  K <- dim(X)[3]
  
  SST <- trace1(Y%*%t(Y))/N/TT
  
  beta <- st                  #starting value for beta-minimization;
  beta_old <- st + 1e5
  obj <- 1e4
  diff_obj <- -1e4
  iter <- 1
  
  XX <- matrix(NA,N*TT,K)
  
  while (max(abs(beta-beta_old))>precision_beta ){
    
    #two stopping criteria for iteration of STEP 1 and STEP 2 below:
    #(1) stop if each component of "beta-beta_old" is smaller than "precision_beta"
    #    ==> this is the good case when we have found a proper local minimum
    #(2) stop if diff_obj>SST*10^-8, i.e. if we made no progress in the objective
    #    function during the last iteration
    #    ==> this is the bad case, where the iteration with that particular
    #    starting value is likely not to converge
  
    #---- STEP 1: CALCULATE FACTORS AND FACTOR LOADINGS FOR GIVEN beta BY PRINCIPAL COMPONENTS: ----
    res <- get_residuals(Y,X,beta)
    temp <- principal_components(res,R)
    lambda <- temp$lambda
    f <- temp$vectors
    res <- res-lambda%*%t(f)                 #residuals after subtracting lambda*f'
    obj_old <- obj                           #save old objective function
    obj <- trace1(t(res)%*%res)/N/T           #LS objective function
    diff_obj <- obj - obj_old                #we hopefully make progress in minimizing objective function, i.e. "diff_obj" should better be negative
    iter <- iter + 1
    
    if(diff_obj<=0){
      
      #---- STEP 2: CALCULATE OLS ESTIMATOR FROM REGRESSING M_lambda*Y*M_f (or just Y) on M_lambda*X*M_f: ----  
      YY <- c(Y)                    #flatten Y, i.e. YY now NTx1 vector
      
      #alternatively, we could define YY as follows, but it should not matter:
      #YY=Y-lambda*(lambda\Y);        #project lambda out of Y
      #YY=( YY'-f*(f\YY') )';         #project f out of Y
      #YY=YY(:);                      #flatten Y, i.e. YY now NTx1 vector       
      
      for (k in 1:K) {
        
        xx <- drop(X[,,k])
        xx <- xx-lambda%*%(qr.solve(lambda,xx))                   #project lambda out of X
        xx <- t(t(xx) - f%*%(qr.solve(f,t(xx))))                  #project f out of X
        XX[,k] <- xx                                             #flatten X, i.e. XX becomes NTxK matrix
        
      }
      
      beta_old <- beta                                    #store old beta
      beta <- solve(t(XX)%*%XX)%*%t(XX)%*%YY               #calculate OLS estimator
      #Here, we use the pseudo-inverse, in case XX"*XX is not invertible to avoide error messages at some points of the optimization.
      #However, at the optimium we should have that beta=XX\YY.

    }
  }
  
  if(diff_obj<=0){
    ef <- 1                 #good solution found
  } else {
    ef <- -1                #no good solution found
  }
  
  obj <- LS_obj(beta,Y,X,R)         #calculate objective function for this beta
  
  return(list(beta=beta,obj=obj,ef=ef,iter=iter))
}


minimize_obj_method2 <- function(Y,X,R,st,precision_beta){
  
  #inputs and outputs as in "minimize_obj_method1"
  temp <- optimize(LS_obj, interval = c(-1e4,1e4), Y=Y, X=X, R=R)
  beta <- temp$minimum
  obj <- temp$objective
  
  return(list(beta=beta,obj=obj))
}


minimize_obj_method3 <- function(Y,X,R,st,precision_beta){
  #inputs and outputs as in "minimize_obj_method1"
  N <- dim(Y)[1]
  TT <- dim(Y)[2]
  K <- dim(X)[3]
  XX <- matrix(NA,N*TT,K)
  
  for(k in 1:K){
    
    xx <- drop(X[,,k])
    XX[,k] <- c(xx)             #flatten X, i.e. XX becomes NTxK matrix (needed for step 2 below)
    
  }
  
  beta <- st                    #starting value for beta-minimization
  beta_old <- st + 1e5
  obj <- 0
  diff_obj <- -1e4 
  iter <- 1
  SST <- trace1(Y%*%t(Y))/N/TT
  obj_save <- 1e4*matrix(1,1,1000)      #we save the objective functions from the previous 1000 iterations
  
  while (max(abs(beta-beta_old))>precision_beta){
    #In case the method does not converge (i.e. if |beta-beta_old| does not become sufficiently small)
    #we need a second stopping criterion, which here we choose relatively conservately, namely, we stop
    #if the objectve function did not improve over the last 1000 iterations by at least 10^-7*trace(Y*Y')/N/T.
                                                           
    #---- STEP 1: CALCULATE FACTORS AND FACTOR LOADINGS FOR GIVEN beta BY PRINCIPAL COMPONENTS: ----  
    res <- get_residuals(Y,X,beta)
    temp <- principal_components(res,R)
    lambda <- temp$lambda
    f <- temp$vectors
    res <- res-lambda%*%t(f)                     #residuals after subtracting lambda*f'
    obj <- trace1(t(res)%*%res)/N/TT                #LS objective function   
    obj_save[(i+1)%%1000+1] <- obj
    diff_obj <- obj - obj_save[(i+1)%%1000+1] 
        #Difference between current objective fct. and objectice fct. from 1000 iterations ago.
        #In this method (as opposed to method 1) this difference is negative by construction.
    
    #---- STEP 2: CALCULATE OLS ESTIMATOR FROM REGRESSING Y-lambda*f' on X: ----  
    YY <- Y-lambda%*%t(f)                #redisuals after proncipal components are subtracted from Y
    YY <- c(YY)                             #flatten Y, i.e. YY now NTx1 vector   
    beta_old <- beta                     #store old beta
    beta <- qr.solve(XX,YY)               #calculate OLS estimator
    
    iter <- iter + 1               #count the number of iterations
  }
  
  if(max(abs(beta-beta_old))<=precision_beta){
    ef <- 1       #good solution found
  } else {
    ef <- -1      #no good solution found
  }
    
  obj <- LS_obj(beta,Y,X,R)        #calculate objective function for this beta
  
  return(list(beta=beta,obj=obj,ef=ef,iter=iter))
}


LS_obj <- function(beta,Y,X,R){
  
  #Calculate LS objective, i.e. SSR/N/T of Y-beta*X after subtracting R largest principal components.
  #INPUT: Kx1 beta, NxT Y, KxNxT X, integer R>=0
  #OUTPUT: scalar objective function
  #
  #COMMENT: within "LS_factor" it is guaranteed that T<=N, so below we diagonalize a TxT matrix
  #(not an NxN) matrix. When using this function outside "LS_factor" one should check whether T<N or N>T 
  #and switch dimensions accordingly, if neccessary.

  res <- get_residuals(Y,X,beta)   
  ev <- eigen(t(res)%*%res)$values
  obj <- sum(ev[1:(dim(Y)[2]-R)])/dim(Y)[1]/dim(Y)[2]
  
  return(obj)
}

get_residuals <- function(Y,X,beta){
  
  #calculate residuals Y-beta*X
  #INPUT: Y = NxT, X = KxNxT, beta=Kx1
  #OUTPUT: res = NxT
  
  res <- Y
  
  for(k in 1:(dim(X)[3])){
    res <- res-beta[k]*drop(X[,,k])
  }
  
  return(res)
}

principal_components <- function(res,R){
  
  #Extract the "R" leading principal components out of the NxT matrix "res".
  #Output: NxR matrix lambda, TxR matrix f (=factor loadings and factors) 
  #
  #COMMENT: within "LS_factor" it is guaranteed that T<=N, so below we diagonalize a TxT matrix (not an NxN) matrix. 
  #When using this function outside "LS_factor" one should check whether T<N or N>T and switch dimensions accordingly, if neccessary.
  
  TT <- dim(res)[2]                       #time dimensions
  eig <- eigen(t(res)%*%res)             #calculate eigenvalues and eigenvectors of TxT matrix
  values <- eig$values
  vectors <- eig$vectors
  
  for(r in 1:R){
    
    vectors[,r] <- vectors[,r]/norm(as.matrix(vectors[,r]),type = "F")          #normalize each principal component (= each factor)
    
    if(mean(vectors[,r]) < 0){
      vectors[,r] <- -vectors[,r]                                     #f(1,r) is normalized to be positive
    }
  }
  
  vectors <- vectors[,1:R]
  lambda <- res%*%vectors              #principal components in the other dimensions (= factor loadings)
  
  return(list(lambda=lambda, vectors=vectors))
}   

trunc1 <- function(A,M1,M2){
  #truncates a symmetric matrix A to the (M1+M2-1) first diagonals
  NN <- dim(A)[1]
  AT <- matrix(0,NN,NN)
  
  for(i in 1:NN){
    for(j in max(i-M1+1,1):min(i+M2-1,NN)){
      AT[i,j] <- A[i,j]
    }
  }
  
  return(AT)
}

trace1 <- function(A) {
  n <- dim(A)[1] # get dimension of matrix
  tr <- 0 # initialize trace value
  
  # Loop over the diagonal elements of the supplied matrix and add the element to tr
  for (k in 1:n) {
    l <- A[k,k]
    tr <- tr + l
  }
  return(tr[[1]])
}





