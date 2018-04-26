rm(list = ls())

library("reshape")
library("car")
library("grid")
library("gridExtra")
library("xtable")
library("plm")
library("ggplot2")

source("LS factor.R")

set.seed(1030)

iter <- 1000
rho0 <- c(0.3,0.9)
T0 <- c(5,10,20,40,80)
N0 <- 100
rhof <- 0.5
sigmaf <- 0.5
table1 <- c()

for(m in 1:length(rho0)){
  temp3 <- c()
  for(n in 1:length(T0)){
    temp1 <- c()
    for(k in 1:iter){
      
      Time <- T0[n]
      e <- matrix(NA,N0,1000+Time)
      lambda <- matrix(NA,N0,1)
      u <- matrix(0,1,1000+Time)
      f <- matrix(0,1,1001+Time)
      Y <- matrix(0,N0,1001+Time)
      
      # Simulation
      
      for(i in 1:N0){
        
        lambda[i] <- rnorm(1,1,1)
        
        for(t in 1:(1000+Time)){
          
          e[i,t] <- rt(n=1,df=5)
          u[t] <- rnorm(1,0,(1-rhof^2)*sigmaf^2)
          f[t+1] <- rhof*f[t]+u[t] 
          Y[i,t+1] <- rho0[m]*Y[i,t]+lambda[i]*f[t+1]+e[i,t]
          
        }
      }
      
      f <- f[,1000:(1000+Time-1)] 
      Y <- as.data.frame(cbind(c(1:nrow(Y)),Y[,1000:(1000+Time-1)])) 
      Y.long <- reshape(Y, varying = list(names(Y)[2:ncol(Y)]),
                           idvar = "V1", timevar = "T", times = c(1:Time), direction = "long")
      
      Ylag <- Y.long[which(Y.long$T != Time),]
      Ylead <- Y.long[which(Y.long$T != 1),]
      
      
      data <- cbind(Ylead, Ylag[,3])
      colnames(data) <- c("id","time","Y","Ylag")
      
      #Estimation
      R <- 2       #specify the number of factors
      
      #OLS
      est.ols <- summary(lm(Y ~ Ylag, data = data))
      rho.ols <- est.ols$coefficients[2,1]
      
      #Dynamic Panel model
      #est.dynamic <- summary(plm(Y ~ Ylag, data = data))
      #rho.dynamic <- est.dynamic$coefficients[1]
      
      #FLS
      Ylag.3d <- as.matrix(Y[,-c(1,ncol(Y))])
      dim(Ylag.3d) <- c(N0,Time-1,1)
      
      est.fls <- LS.factor(as.matrix(Y[,-c(1:2)]), Ylag.3d, R, method="m1", repMIN=5, repMAX=5, report=FALSE)
      rho.FLS <- est.fls$bcorr3 + rho0[m]
      
      #BC_FLS
      Ylag.3d <- as.matrix(Y[,-c(1,ncol(Y))])
      dim(Ylag.3d) <- c(N0,Time-1,1)
      
      est.fls <- LS.factor(as.matrix(Y[,-c(1:2)]), Ylag.3d, R, method="m1", repMIN=5, repMAX=5, report=FALSE)
      rho.BC_FLS <- est.fls$beta
      
      temp1 <- rbind(temp1,c(rho.ols, rho.FLS, rho.BC_FLS))
      mean_bias <- colMeans(temp1) - rho0[m]
      std <- apply(temp1,2,sd)
      rmse <- sqrt(colMeans((temp1 - rho0[m])^2))
      temp2 <- rbind(mean_bias,std,rmse)
    }
    temp3 <- rbind(temp3,temp2)
    cat("Completed rho0 = ", rho0[m], "T= ", Time, "\n")
  }
  table1 <- cbind(table1,temp3)
}

colnames(table1) <- c("OLS","FLS","BC-FLS","OLS","FLS","BC-FLS")
# Create PDF file
pdf("Table 2.pdf", width = 12, height = 12, onefile = TRUE)
grid.table(format(round(table1,digits = 4), nsmall = 4))
dev.off()   
  
# Create Latex code
xtable(table1)
