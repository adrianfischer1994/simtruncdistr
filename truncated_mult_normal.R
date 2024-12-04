
library(parallel)
numCores <- 10
cl <- makeCluster(numCores)

clusterEvalQ(cl,{
  
  library(MASS)
  library(cubature)
  library(matrixcalc)
  
  set.seed(1)
  n <- 100
  
  #calculation of the normalizing constant for cuboids
  norm_c <- function(mu,sigmainv,a) {
    d <- length(a[1,])
    func_int <- function(x) {
      return(exp(-0.5*t(x-mu)%*%sigmainv%*%(x-mu)))
    }
    return(cubintegrate(func_int,lower=a[1,],upper=a[2,])$integral)
  }
  
})

sim <- function(param) {
  
  
  mu <- param[[1]]
  sigma <- param[[2]]
  d <- length(mu)
  
  
  #functions for Stein estimation ellipse
  # a <- c(1,4,1,2)
  # kappa <- function(x) {
  #   return(sum(x^2/a^2)-1)
  # }
  # dkappa <- function(x) {
  #   return(2*x/a^2)
  # }
  
  #functions for Stein estimation cuboid
  a <- rbind(rep(-1,d),rep(1,d))
  kappa <- function(x) {
    return(prod((x-a[1,])*(x-a[2,])))
  }
  dkappa <- function(x) {
    dim <- length(x)
    erg <- rep(NA,dim)
    for(i in 1:dim) {
      a_t <- as.matrix(a[,-i])
      x_t <- x[-i]
      erg[i] <- (2*x[i]-sum(a[,i]))*prod((x_t-a_t[1,])*(x_t-a_t[2,]))
    }
    return(erg)
  }
  
  MLE <- list()
  SM <- list()
  ST <- list()
  
  MLEtime <- rep(NA,m)
  SMtime <- rep(NA,m)
  STtime <- rep(NA,m)
  
  assign("MLE_nonex", 0, env=globalenv())
  assign("SM_nonex", 0, env=globalenv())
  assign("ST_nonex", 0, env=globalenv())
  
  for(i in 1:m) {
    
    X <- matrix(NA,n,d)
    n_tmp <- 1
    
    #Sampling for truncation wrt ellipse
    # while(n_tmp<(n+1)) {
    #   x_tmp <- mvrnorm(1,mu,sigma)
    #   if(sum(x_tmp^2/a^2)<1) {
    #     X[n_tmp,] <- x_tmp
    #     n_tmp <- n_tmp+1
    #   }
    # }
    
    #Sampling for truncation wrt cuboid
    while(n_tmp<(n+1)) {
      x_tmp <- mvrnorm(1,mu,sigma)
      in_dom <- TRUE
      for(j in 1:d) {
        if(x_tmp[j]<a[1,j] || x_tmp[j]>a[2,j]) {
          in_dom <- FALSE
        }
      }
      if(in_dom) {
        X[n_tmp,] <- x_tmp
        n_tmp <- n_tmp+1
      }
    }
  
    
    #MLE
    tryCatch({
      starttime <- Sys.time()
      ll <- function(par) {
        d <- length(a[1,])
        mu <- par[1:d]
        sigma <- matrix(0,d,d)
        for(j in 1:d) {
          for(k in 1:j) {
            sigma[j,k] <- par[d+((j-1)^2+j-1)/2+k]
          }
        }
        sigma <- sigma%*%t(sigma)
        sigmainv <- solve(sigma)
        C <- norm_c(mu,sigmainv,a)
        func_temp <- function(x) t(x-mu)%*%sigmainv%*%(x-mu)
        return(n*log(C)+0.5*sum(apply(X,1,func_temp)))
      }
      MLEstart <- numeric(d+d*(d+1)/2)
      for(j in 1:d) {
        MLEstart[d+((j-1)^2+j-1)/2+j] <- 1
      }
      MLEtemp <- R.utils::withTimeout( optim(MLEstart,ll)$par , timeout = 20, onTimeout = "silent")
      MLEmu <- MLEtemp[1:d]
      MLEsigma <- matrix(0,d,d)
      for(j in 1:d) {
        for(k in 1:j) {
          MLEsigma[j,k] <- MLEtemp[d+((j-1)^2+j-1)/2+k]
        }
      }
      MLEsigma <- MLEsigma%*%t(MLEsigma)
      if(all(!is.nan(MLEmu)) & all(!is.na(MLEmu)) & all(!is.nan(MLEsigma)) & all(!is.na(MLEsigma)) & all(eigen(MLEsigma)$values>0)) {
        MLE[[length(MLE)+1]] <- list(MLEmu,MLEsigma)
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
    
    
    #SM ##does only work for cuboids
    tryCatch({
      starttime <- Sys.time()
      smf <- function(par) {
        d <- length(a[1,])
        mu <- par[1:d]
        sigma <- matrix(0,d,d)
        for(j in 1:d) {
          for(k in 1:j) {
            sigma[j,k] <- par[d+((j-1)^2+j-1)/2+k]
          }
        }
        sigma <- sigma%*%t(sigma)
        sigmainv <- solve(sigma)
        func_temp <- function(x) sum(pmin(x-a[1,],a[2,]-x)*(sigmainv%*%(x-mu))^2)
        func_temp2 <- function(x) sum((-1*(x<(a[1,]+a[2,])/2)+1*(x>(a[1,]+a[2,])/2))*sigmainv%*%(x-mu)-pmin(x-a[1,],a[2,]-x)*diag(sigmainv))
        return(mean(apply(X,1,func_temp))+2*mean(apply(X,1,func_temp2)))
      }
      SMstart <- numeric(d+d*(d+1)/2)
      for(j in 1:d) {
        SMstart[d+((j-1)^2+j-1)/2+j] <- 1
      }
      SMtemp <-  R.utils::withTimeout( optim(SMstart,smf)$par , timeout = 20, onTimeout = "silent")
      SMmu <- SMtemp[1:d]
      SMsigma <- matrix(0,d,d)
      for(j in 1:d) {
        for(k in 1:j) {
          SMsigma[j,k] <- SMtemp[d+((j-1)^2+j-1)/2+k]
        }
      }
      SMsigma <- SMsigma%*%t(SMsigma)
      if(all(!is.nan(SMmu)) & all(!is.na(SMmu)) & all(!is.nan(SMsigma)) & all(!is.na(SMsigma)) & all(eigen(SMsigma)$values>0)) {
        SM[[length(SM)+1]] <- list(SMmu,SMsigma)
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
      f1 <- mean(apply(X,1,kappa))

      xf1 <- rep(0,d)
      for(j in 1:d) {
        xf1[j] <- mean(X[,j]*apply(X,1,kappa))
      }

      dkp_tmp <- apply(X,1,dkappa)
      df1 <- rowMeans(dkp_tmp)

      f2 <- xf1

      xf2 <- matrix(0,d,d)
      for(j in 1:d) {
        for(k in 1:d) {
          xf2[j,k] <- mean(X[,j]*X[,k]*apply(X,1,kappa))
        }
      }

      df2 <- matrix(0,d,d)
      for(j in 1:d) {
        for(k in 1:d) {
          if(j==k) {
            df2[j,k] <- mean(apply(X,1,kappa))+mean(X[,j]*dkp_tmp[j,])
          } else {
            df2[j,k] <- mean(X[,k]*dkp_tmp[j,])
          }
        }
      }


      STsigma <- (t(xf2)*f1-xf1%*%t(f2))%*%solve(df2*f1-df1%*%t(f2))
      STsigma <- 0.5*(t(STsigma)+STsigma)
      STmu <- (xf1-STsigma%*%df1)/f1
      
      if(all(!is.nan(STmu)) & all(!is.na(STmu)) & all(!is.nan(STsigma)) & all(!is.na(STsigma)) & all(eigen(STsigma)$values>0)) {
        ST[[length(ST)+1]] <- list(STmu,STsigma)
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
  
  mse <- matrix(0,3,2)
  ex <- numeric(3)
  time <- numeric(3)

  l1 <- length(MLE)
  l2 <- length(SM)
  l3 <- length(ST)

  if(l1>0) {
    MLEmu_n <- numeric(l1)
    MLEsigma_n <- numeric(l1)
    
    for(i in 1:l1) {
      MLEmu_n[i] <- norm(MLE[[i]][[1]]-mu, type="2")
      MLEsigma_n[i] <- spectral.norm(MLE[[i]][[2]]-sigma)
    }
    
    mse[1,] <- c(mean(MLEmu_n),mean(MLEsigma_n))
    time[1] <- mean(MLEtime, na.rm=TRUE)
  } else {
    mse[1,] <- NaN
    time[1] <- NaN
  }
  ex[1] <- MLE_nonex
  
  if(l2>0) {
    SMmu_n <- numeric(l2)
    SMsigma_n <- numeric(l2)
    
    for(i in 1:l2) {
      SMmu_n[i] <- norm(SM[[i]][[1]]-mu, type="2")
      SMsigma_n[i] <- spectral.norm(SM[[i]][[2]]-sigma)
    }
    
    mse[2,] <- c(mean(SMmu_n),mean(SMsigma_n))
    time[2] <- mean(SMtime, na.rm=TRUE)
  } else {
    mse[2,] <- NaN
    time[2] <- NaN
  }
  ex[2] <- SM_nonex
  
  if(l3>0) {
    STmu_n <- numeric(l3)
    STsigma_n <- numeric(l3)
    
    for(i in 1:l3) {
      STmu_n[i] <- norm(ST[[i]][[1]]-mu, type="2")
      STsigma_n[i] <- spectral.norm(ST[[i]][[2]]-sigma)
    }
    
    mse[3,] <- c(mean(STmu_n),mean(STsigma_n))
    time[3] <- mean(STtime, na.rm=TRUE)
  } else {
    mse[3,] <- NaN
    time[3] <- NaN
  }
  ex[3] <- ST_nonex
  
  return(list(mse,ex,time))
}

m <- 2
clusterExport(cl,c("m"))

mu1 <- c(0,0)
sigma1 <- diag(2)

mu2 <- c(0.5,0.5)
sigma2 <- diag(2)

mu3 <- c(0.5,0.5)
sigma3 <- 0.5*diag(2)

mu4 <- c(0,0)
sigma4 <- 2*diag(2)

mu5 <- c(0,0)
sigma5 <- 0.2*diag(2)

mu6 <- c(0.8,-0.2)
sigma6 <- 0.5*diag(2)

mu7 <- c(0,0)
sigma7 <- cbind(c(0.5,0.4),c(0.4,0.5))

mu8 <- c(0,0)
sigma8 <- cbind(c(0.8,-0.7),c(-0.7,0.9))

mu9 <- c(0.3,-0.2)
sigma9 <- cbind(c(0.2,0.1),c(0.1,0.4))

mu10 <- c(0.5,0.5)
sigma10 <- cbind(c(0.1,0.1),c(0.1,0.8))

param <- list(list(mu1,sigma1),list(mu2,sigma2),list(mu3,sigma3),list(mu4,sigma4),
              list(mu5,sigma5),list(mu6,sigma6),list(mu7,sigma7),list(mu8,sigma8),list(mu9,sigma9),list(mu10,sigma10))

erg <- parLapply(cl,param,sim)

stopCluster(cl)

saveRDS(erg, file="result.RData")

nb_est <- length(erg[[1]][[1]][,1]) #number of estimators to compare
nb_par <- length(erg[[1]][[1]][1,]) #number of unknown parameters
nb_parco <- length(erg) #number of parameter constellations

trans <- matrix(0,length(erg),nb_est*nb_par)
for(i in 1:length(erg)) {
  for(j in 1:nb_par) {
    for(k in 1:nb_est) {
      trans[i,(j-1)*nb_est+k] <- erg[[i]][[1]][k,j]
    }    
  }
}
trans_nonex <- matrix(0,nb_parco,nb_est)
for(i in 1:nb_parco) {
  for(j in 1:nb_est) {
    trans_nonex[i,j] <- erg[[i]][[2]][j]
  }
}
trans_time <- matrix(0,nb_parco,nb_est)
for(i in 1:nb_parco) {
  for(j in 1:nb_est) {
    trans_time[i,j] <- erg[[i]][[3]][j]
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
for(i in 1:(nb_est*1+2)) {
  table <- paste(table,"c",sep="")
}
table <- paste(table,"}\n $\\mu_0,\\Sigma_0$ & & \\multicolumn{",nb_est,"}{c}{MSE} \\\\ \n",sep="")
for(i in 1:length(erg)) {
  d <- length(param[[i]][[1]])
  pt1 <- ""
  for(j in 1:d) {
    if(j==d) {
      pt1 <- paste(pt1,param[[i]][[1]][j],sep="")
    } else {
      pt1 <- paste(pt1,param[[i]][[1]][j]," \\\\ ",sep="")
    }
  }
  pt2 <- ""
  for(j in 1:d) {
    for(k in 1:d) {
      if(k==d) {
        if(j==d) {
          pt2 <- paste(pt2,param[[i]][[2]][j,k],sep="")
        } else {
          pt2 <- paste(pt2,param[[i]][[2]][j,k]," \\\\ ",sep="")
        }
      } else {
        pt2 <- paste(pt2,param[[i]][[2]][j,k]," & ",sep="")
      }
    }
  }
  table <- paste(table,"\\multirow{2}{*}{$\\begin{pmatrix}",pt1,"\\end{pmatrix},\\begin{pmatrix}",pt2,"\\end{pmatrix}$} & $\\mu$\\Tstrut ",sep="")
  for(j in 1:length(trans[1,])) { 
    ent <- rnd(trans[i,j])
    if(j %% (nb_est) == 0) {
      if(j == nb_est){
        table <- paste(table," & $",ent,"$", sep="")
        table <- paste(table," \\\\ & $\\Sigma$\\Bstrut ", sep="")
      } else {
        table <- paste(table," & $",ent,"$ \\\\", sep="")
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
table <- paste(table,"}\n $\\mu_0,\\Sigma_0$ & & \\multicolumn{",nb_est,"}{c}{NE} & \\multicolumn{",nb_est,"}{c}{Time} \\\\ \n",sep="")
for(i in 1:length(erg)) {
  d <- length(param[[i]][[1]])
  pt1 <- ""
  for(j in 1:d) {
    if(j==d) {
      pt1 <- paste(pt1,param[[i]][[1]][j],sep="")
    } else {
      pt1 <- paste(pt1,param[[i]][[1]][j]," \\\\ ",sep="")
    }
  }
  pt2 <- ""
  for(j in 1:d) {
    for(k in 1:d) {
      if(k==d) {
        if(j==d) {
          pt2 <- paste(pt2,param[[i]][[2]][j,k],sep="")
        } else {
          pt2 <- paste(pt2,param[[i]][[2]][j,k]," \\\\ ",sep="")
        }
      } else {
        pt2 <- paste(pt2,param[[i]][[2]][j,k]," & ",sep="")
      }
    }
  }
  table <- paste(table,"\\multirow{2}{*}{$\\begin{pmatrix}",pt1,"\\end{pmatrix},\\begin{pmatrix}",pt2,"\\end{pmatrix}$} & $\\mu$\\Tstrut ",sep="")
  for(k in 1:nb_est) {
    table <- paste(table," & \\multirow{",nb_par,"}{*}{$",round(trans_nonex[i,k]/m,2)*100,"$}", sep="")
  }
  for(k in 1:nb_est) {
    table <- paste(table," & \\multirow{",nb_par,"}{*}{$",rnd(trans_time[i,k]),"$}", sep="")
  }
  table <- paste(table," \\\\ & $\\Sigma$\\Bstrut & & & & & & \\\\ ", sep="")
  table <- paste(table," \\hline \n", sep="")
}

table <- paste(table,"\\end{tabular} \n\\end{table}", sep="")

sink("table2.txt")
cat(table)
sink()
