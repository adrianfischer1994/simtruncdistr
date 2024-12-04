library(parallel)

numCores <- 7
cl <- makeCluster(numCores)

clusterEvalQ(cl,{
  
  set.seed(1)
  n <- 500
  
  #calculation of the normalizing constant for circles
  norm_c <- function(mu,sigma,alpha,beta,a,r) {
    
    inner_int <- function(y) {
      inner_func <- function(x) {
        return(dnorm(x,mean=mu,sd=sqrt(sigma))*dbeta(y,shape1=alpha,shape2=beta))
      }
      inner_func <- Vectorize(inner_func)
      return(integrate(inner_func,lower=-sqrt(r^2-(y-a[2])^2)+a[1],upper=sqrt(r^2-(y-a[2])^2)+a[1])$value)
    }
    inner_int <- Vectorize(inner_int)
    
    return(integrate(inner_int,lower=a[2]-r,upper=a[2]+r)$value)
  }
  
})

sim <- function(param) {
  
  #sigma corresponds to sigma squared here
  mu <- param[1]
  sigma <- param[2]
  alpha <- param[3]
  beta <- param[4]
  
  #functions for Stein estimation for circles
  r <- 0.5
  a <- c(0,0.5)
  kappa <- function(x) {
    return((x[1]-a[1])^2+(x[2]-a[2])^2-r^2)
  }
  dkappa <- function(x) {
    return(2*(x-a))
  }
  
  f1 <- function(x) {
    return(kappa(x))
  }
  df1 <- function(x) {
    return(dkappa(x))
  }
  f2 <- function(x) {
    return((x[1]+x[2])*kappa(x))
  }
  df2 <- function(x) {
    return(c(1,1)*kappa(x)+(x[1]+x[2])*dkappa(x))
  }
  
  
  MLE <- matrix(NA,m,4)
  SM <- matrix(NA,m,4)
  ST <- matrix(NA,m,4)
  
  MLEtime <- rep(NA,m)
  SMtime <- rep(NA,m)
  STtime <- rep(NA,m)
  
  assign("MLE_nonex", 0, env=globalenv())
  assign("SM_nonex", 0, env=globalenv())
  assign("ST_nonex", 0, env=globalenv())
  
  for(i in 1:m) {
    
    X <- matrix(NA,n,2)
    n_tmp <- 1
    
    #Sampling for truncation wrt circles
    while(n_tmp<(n+1)) {
      x_tmp <- c(rnorm(1,mu,sqrt(sigma)),rbeta(1,shape1=alpha,shape2=beta))
      if((x_tmp[1]-a[1])^2+(x_tmp[2]-a[2])^2<r^2) {
        X[n_tmp,] <- x_tmp
        n_tmp <- n_tmp+1
      }
    }
              
    #MLE
    tryCatch({
      starttime <- Sys.time()
      ll <- function(par) {
        C <- norm_c(par[1],par[2],par[3],par[4],a,r)
        return(n*log(C)-sum(log(dnorm(X[,1],mean=par[1],sd=sqrt(par[2]))))-sum(log(dbeta(X[,2],shape1=par[3],shape2=par[4]))))
      }

      MLEparam <- R.utils::withTimeout( optim(c(0,1,1,1),ll,method="L-BFGS-B",lower = c(-Inf,0,0,0), upper = c(Inf,Inf,Inf,Inf))$par , timeout = 20, onTimeout = "silent")

      if(all(!is.nan(MLEparam)) & all(!is.na(MLEparam)) & all(!is.null(MLEparam)) & MLEparam[2]>0 & MLEparam[3]>0 & MLEparam[4]>0) {
        MLE[i,] <- MLEparam
        endtime <- Sys.time()
        MLEtime[i] <- as.double(endtime-starttime, units = "secs")
      } else {
        nonex <- get("MLE_nonex", env=globalenv())
        assign("MLE_nonex", nonex+1, env=globalenv())
      }
    },error=function(cond) {
      nonex <- get("MLE_nonex", env=globalenv())
      assign("MLE_nonex", nonex+1, env=globalenv())
    })
    
    #SM
    tryCatch({
      starttime <- Sys.time()
      
      #This is the numericl optimisation code, but explicit estimators are available
      # smf <- function(par) {
      #   func_temp <- function(x) (r-norm(x-a,type="2"))*((x[1]-par[1])^2/par[2]^2+((par[3]-1)/x[2]-(par[4]-1)/(1-x[2]))^2) + 2*(r-norm(x-a,type="2"))*(-1/par[2]-(par[3]-1)/x[2]^2-(par[4]-1)/(1-x[2])^2) + 
      #     2*((x[1]-a[1])*((x[1]-par[1])/par[2])-(x[2]-a[2])*((par[3]-1)/x[2]-(par[4]-1)/(1-x[2])))/norm(x-a,type="2")
      #   return(mean(apply(X,1,func_temp)))
      # }
      # SMparam <- optim(c(0,1,1,1),smf)$par 
      
      func_temp <- function(x) (r-norm(x-a,type="2"))
      M1 <- mean(apply(X,1,func_temp))
      
      func_temp <- function(x) -2*(r-norm(x-a,type="2"))*x[1]
      M2 <- mean(apply(X,1,func_temp))
      
      func_temp <- function(x) (r-norm(x-a,type="2"))*x[1]^2
      M3 <- mean(apply(X,1,func_temp))
      
      func_temp <- function(x) (r-norm(x-a,type="2"))/x[2]^2
      M4 <- mean(apply(X,1,func_temp))
      
      func_temp <- function(x) (r-norm(x-a,type="2"))/(1-x[2])^2
      M5 <- mean(apply(X,1,func_temp))
      
      func_temp <- function(x) -(r-norm(x-a,type="2"))*2/(x[2]*(1-x[2]))
      M6 <- mean(apply(X,1,func_temp))
      
      func_temp <- function(x) -2*(r-norm(x-a,type="2"))+2*(x[1]-a[1])/norm(x-a,type="2")*x[1]
      M7 <- mean(apply(X,1,func_temp))
      
      func_temp <- function(x) -2*(x[1]-a[1])/norm(x-a,type="2")
      M8 <- mean(apply(X,1,func_temp))
      
      func_temp <- function(x) (r-norm(x-a,type="2"))*2*1/x[2]*(1/(1-x[2])-1/x[2]) - 2*(r-norm(x-a,type="2"))*1/x[2]^2 - 2*(x[2]-a[2])/norm(x-a,type="2")*1/x[2]
      M9 <- mean(apply(X,1,func_temp))
      
      func_temp <- function(x) -(r-norm(x-a,type="2"))*2*1/(1-x[2])*(1/(1-x[2])-1/x[2]) - 2*(r-norm(x-a,type="2"))*1/(1-x[2])^2 + 2*(x[2]-a[2])/norm(x-a,type="2")*1/(1-x[2])
      M10 <- mean(apply(X,1,func_temp))
      
      SMparam <- rep(NA,4)
      SMparam[1] <- (2*M3*M8-M2*M7)/(2*M1*M7-M2*M8)
      SMparam[2] <- (M2^2-4*M1*M3)/(2*M1*M7-M2*M8)
      SMparam[3] <- (M6*M10-2*M5*M9)/(4*M4*M5-M6^2)
      SMparam[4] <- (M6*M9-2*M4*M10)/(4*M4*M5-M6^2)
      
      if(all(!is.nan(SMparam)) & all(!is.na(SMparam)) & all(!is.null(SMparam)) & SMparam[2]>0 & SMparam[3]>0 & SMparam[4]>0) {
        SM[i,] <- SMparam
        endtime <- Sys.time()
        SMtime[i] <- as.double(endtime-starttime, units = "secs")
      } else {
        nonex <- get("SM_nonex", env=globalenv())
        assign("SM_nonex", nonex+1, env=globalenv())
      }
    },error=function(cond) {
      nonex <- get("SM_nonex", env=globalenv())
      assign("SM_nonex", nonex+1, env=globalenv())
    })
    
    #ST
    tryCatch({
      starttime <- Sys.time()
      
      f1x <- apply(X,1,f1)
      df1x <- apply(X,1,df1)
      f2x <- apply(X,1,f2)
      df2x <- apply(X,1,df2)

      M1 <- mean(X[,2]*(1-X[,2])*df2x[1,])
      M2 <- mean(X[,2]*(1-X[,2])*X[,1]*f1x)
      M3 <- mean(X[,2]*(1-X[,2])*df1x[1,])
      M4 <- mean(X[,2]*(1-X[,2])*X[,1]*f2x)
      M5 <- mean(X[,2]*(1-X[,2])*f1x)
      M6 <- mean(X[,2]*(1-X[,2])*f2x)

      STmu <- (M1*M2-M3*M4)/(M5*M1-M3*M6)
      STsigma <- (M5*M4-M6*M2)/(M5*M1-M3*M6)

      N1 <- mean(X[,2]*f1x)
      N2 <- mean(X[,2]*f2x)
      N3 <- mean((1-X[,2])*X[,2]*df1x[2,])
      N4 <- mean((1-X[,2])*X[,2]*df2x[2,])
      N5 <- mean((X[,2]-1)*f1x)
      N6 <- mean((X[,2]-1)*f2x)

      STalpha <- (N1*N4-N2*N3)/(N1*N6-N2*N5)
      STbeta <- (N6*N3-N5*N4)/(N1*N6-N2*N5)

      STparam <- c(STmu,STsigma,STalpha,STbeta)

      if(all(!is.nan(STparam)) & all(!is.na(STparam)) & STsigma>0 & STalpha>0 & STbeta>0) {
        ST[i,] <- STparam
        endtime <- Sys.time()
        STtime[i] <- as.double(endtime-starttime, units = "secs")
      } else {
        nonex <- get("ST_nonex", env=globalenv())
        assign("ST_nonex", nonex+1, env=globalenv())
      }

    },error=function(cond) {
      nonex <- get("ST_nonex", env=globalenv())
      assign("ST_nonex", nonex+1, env=globalenv())
    })
  }
  
  bias <- matrix(0,3,4)
  mse <- matrix(0,3,4)
  ex <- numeric(3)
  time <- numeric(3)
  
  bias[1,] <- c(mean(MLE[,1]-mu,na.rm=TRUE),mean(MLE[,2]-sigma^2,na.rm=TRUE),mean(MLE[,3]-alpha,na.rm=TRUE),mean(MLE[,4]-beta,na.rm=TRUE))
  bias[2,] <- c(mean(SM[,1]-mu,na.rm=TRUE),mean(SM[,2]-sigma,na.rm=TRUE),mean(SM[,3]-alpha,na.rm=TRUE),mean(SM[,4]-beta,na.rm=TRUE))
  bias[3,] <- c(mean(ST[,1]-mu,na.rm=TRUE),mean(ST[,2]-sigma,na.rm=TRUE),mean(ST[,3]-alpha,na.rm=TRUE),mean(ST[,4]-beta,na.rm=TRUE))
  
  mse[1,] <- c(mean((MLE[,1]-mu)^2,na.rm=TRUE),mean((MLE[,2]-sigma^2)^2,na.rm=TRUE),mean((MLE[,3]-alpha)^2,na.rm=TRUE),mean((MLE[,4]-beta)^2,na.rm=TRUE))
  mse[2,] <- c(mean((SM[,1]-mu)^2,na.rm=TRUE),mean((SM[,2]-sigma^2)^2,na.rm=TRUE),mean((SM[,3]-alpha)^2,na.rm=TRUE),mean((SM[,4]-beta)^2,na.rm=TRUE))
  mse[3,] <- c(mean((ST[,1]-mu)^2,na.rm=TRUE),mean((ST[,2]-sigma^2)^2,na.rm=TRUE),mean((ST[,3]-alpha)^2,na.rm=TRUE),mean((ST[,4]-beta)^2,na.rm=TRUE))
  
  ex[1] <- MLE_nonex
  ex[2] <- SM_nonex
  ex[3] <- ST_nonex
  
  time[1] <- mean(MLEtime, na.rm=TRUE)
  time[2] <- mean(SMtime, na.rm=TRUE)
  time[3] <- mean(STtime, na.rm=TRUE)
  
  return(list(bias,mse,ex,time))
}

m <- 3
clusterExport(cl,c("m"))

param <- list(c(1,2,1,1),c(0.5,0.1,4,5),c(0,1,1,1.5),c(0,0.1,0.5,3),c(0.2,0.3,0.1,0.4),c(0,1.5,1,0.5),c(0,0.4,2,2))

erg <- parLapply(cl,param,sim)

stopCluster(cl)

saveRDS(erg, file="result.RData")

nb_est <- length(erg[[1]][[1]][,1]) #number of estimators to compare
nb_par <- length(erg[[1]][[1]][1,]) #number of unknown parameters

trans <- matrix(0,length(erg),nb_est*nb_par*2)
for(i in 1:length(erg)) {
  for(j in 1:nb_par) {
    for(l in 1:2) {
      for(k in 1:nb_est) {
        trans[i,(j-1)*2*nb_est+(l-1)*nb_est+k] <- erg[[i]][[l]][k,j]
      }    
    }
  }
}

trans_nonex <- matrix(0,length(param),nb_est)
for(i in 1:length(param)) {
  for(j in 1:nb_est) {
    trans_nonex[i,j] <- erg[[i]][[3]][j]
  }
}
trans_time <- matrix(0,length(param),nb_est)
for(i in 1:length(param)) {
  for(j in 1:nb_est) {
    trans_time[i,j] <- erg[[i]][[4]][j]
  }
}

rnd <- function(ent) {
  if(is.nan(ent)) {
    return("NaN")
  }
  else { 
    if(abs(ent) < 0.01 | abs(ent) > 10000) {
      ex <- ifelse(ent == 0, 0, floor(log10(abs(ent))))
      mant <- round(ent/10^ex,2)
      ent <- paste(mant,"\\text{e",ex,"}",sep="")
    } else {
      if(1<abs(ent) & abs(ent)<10) {
        ent <- round(ent,2)
      } else {
        if(10<abs(ent) & abs(ent)<100) {
          ent <- round(ent,1)
        } else {
          if(100<abs(ent) & abs(ent)<10000) {
            ent <- round(ent,0)
          } else {
            ent <- round(ent,3)
          }
        }
      }
    }
    return(ent)
  }
}

table <- "\\begin{table} \n\\centering\n\\begin{tabular}{"
for(i in 1:(nb_est*2+2)) {
  table <- paste(table,"c",sep="")
}
table <- paste(table,"}\n $\\theta_0$ & & \\multicolumn{",nb_est,"}{c}{Bias} & \\multicolumn{",nb_est,"}{c}{MSE} \\\\ \n",sep="")
for(i in 1:length(erg)) {
  d <- length(param[[i]])
  pt1 <- ""
  for(j in 1:d) {
    if(j==d) {
      pt1 <- paste(pt1,param[[i]][j],sep="")
    } else {
      pt1 <- paste(pt1,param[[i]][j]," \\\\ ",sep="")
    }
  }
  table <- paste(table,"\\multirow{4}{*}{$\\begin{pmatrix}",pt1,"\\end{pmatrix}$} & $\\mu$\\Tstrut ",sep="")
  for(j in 1:length(trans[1,])) { 
    ent <- rnd(trans[i,j])
    if(j %% (2*nb_est) == 0) {
      if(j == (2*nb_est)){
        table <- paste(table," & $",ent,"$", sep="")
        table <- paste(table," \\\\ & $\\sigma^2$ ", sep="")
      } else {
        if(j == (4*nb_est)) {
          table <- paste(table," & $",ent,"$ \\\\ & $\\alpha$ ", sep="")
        } else {
          if(j == (6*nb_est)) {
            table <- paste(table," & $",ent,"$ \\\\ & $\\beta$\\Bstrut ", sep="")
          } else {
            table <- paste(table," & $",ent,"$ \\\\", sep="")
          }
        }
      }
    } else {
      table <- paste(table," & $",ent,"$", sep="")
    }
  }
  table <- paste(table," \\hline \n", sep="")
}

table <- paste(table,"\\end{tabular} \n\\end{table}", sep="")

sink("table1.txt")
cat(table)
sink()


table <- "\\begin{table} \n\\centering\n\\begin{tabular}{"
for(i in 1:(nb_est*2+2)) {
  table <- paste(table,"c",sep="")
}
table <- paste(table,"}\n $\\theta_0$ & & \\multicolumn{",nb_est,"}{c}{NE} & \\multicolumn{",nb_est,"}{c}{Time} \\\\ \n",sep="")
for(i in 1:length(erg)) {
  d <- length(param[[i]])
  pt1 <- ""
  for(j in 1:d) {
    if(j==d) {
      pt1 <- paste(pt1,param[[i]][j],sep="")
    } else {
      pt1 <- paste(pt1,param[[i]][j]," \\\\ ",sep="")
    }
  }
  table <- paste(table,"\\multirow{4}{*}{$\\begin{pmatrix}",pt1,"\\end{pmatrix}$} & $\\mu$\\Tstrut ",sep="")
  for(k in 1:nb_est) {
    table <- paste(table," & \\multirow{",nb_par,"}{*}{$",round(trans_nonex[i,k]/m,2)*100,"$}", sep="")
  }
  for(k in 1:nb_est) {
    table <- paste(table," & \\multirow{",nb_par,"}{*}{$",rnd(trans_time[i,k]),"$}", sep="")
  }
  table <- paste(table," \\\\ & $\\sigma^2$ & & & & & & \\\\ & $\\alpha$ & & & & & & \\\\ & $\\beta\\Bstrut$ & & & & & & \\\\ ", sep="")
  table <- paste(table," \\hline \n", sep="")
}

table <- paste(table,"\\end{tabular} \n\\end{table}", sep="")

sink("table2.txt")
cat(table)
sink()