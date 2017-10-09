
resampling.test.stat.kins_adj.normal_appx <- function(obj.res,obj.list,Z1,resample){
  
  n.obj <- length(obj.list)
  n.pheno <- obj.res$n.pheno
  n <- obj.res$n.all
  m <- ncol(Z1)
  
  V_y <- obj.res$y.Cov
  V_p <- obj.res$pc.Cov
  L <- obj.res$L

  R <- R1 <- A <- Lambda <- Q <- wt <- vector("list",n.obj)
  for(i in 1:n.obj){
    R[[i]] <- obj.list[[i]]$R
    R1[[i]] <- mat.sqrt(R[[i]])
  }

  D <- obj.res$D
  d <- diag(D)
  V_g <- obj.res$V_g
  
  res.mat = matrix(obj.res$res, byrow=TRUE, ncol=n.pheno)
  res <- obj.res$res
  V.item.inv = obj.res$V.item.inv
  V.item <- solve(V.item.inv)
  
  Vcurl <- list()
  for(i in 1:n){
    Vcurl[[i]] <-mat.sqrt(solve(V.item + d[i]*V_g)) 
  }
  
  U <- obj.res$U
  Z <- Z1
  Z1 <- t(U)%*%Z1
  
  Qs <- matrix(0,ncol = n.obj,nrow = resample)
  
  for(j in 1:resample){
    
    res.1 <- scale(as.matrix(rnorm(n*n.pheno),ncol = n.pheno),center = TRUE,scale = F)
#   res.1 <- as.matrix(rnorm(n*n.pheno),ncol = n.pheno)
   res.1 <- as.vector(t(res.1))
    
    res.Vhalf<-rep(0,n*n.pheno)
    for(id1 in 1:n){
      idx<-((id1-1)*n.pheno+1):(id1*n.pheno)
      
      res.Vhalf[idx]<- as.matrix(Vcurl[[id1]])%*% res.1[idx]
      
    }
    res.mat.2 = scale(matrix(res.Vhalf, byrow=TRUE, ncol=n.pheno),center = TRUE,scale = F)
    #   res.mat.2 <- scale(obj.res$U%*%res.mat.2,center = T,scale = F)
    
    
    Q.temp1<-t(res.mat.2) %*% Z1

   trs <- vector("list",n.obj)
    for(indx in 1:n.obj){
	trs[[indx]] <- matrix(0,ncol = m, nrow = n.pheno)
    }
    
    tr1 <- tr2 <- tr3 <- tr4 <- matrix(0,ncol = m, nrow = n.pheno)
    for(indx in 1:n.obj){
     for(i in 1:m){
 	     trs[[indx]][,i] <- R1[[indx]]%*%Q.temp1[,i];
   	}
    }
    for(indx in 1:n.obj){
	Qs[j,indx] <- sum(trs[[indx]]^2)
    }
    m1 <- m*n.pheno
  }
  
  return(list(Test.Stat = Qs,Sigma_Ps = R,dim.W = m1))
}



resampling.pvalues.kins_adj <- function(obj.res,obj.list,Z1,resample){
  
  n.obj <- length(obj.list);
  n.pheno <- obj.res$n.pheno
  n <- obj.res$n.all
  m <- ncol(Z1)
  
  V_y <- obj.res$y.Cov
  V_p <- obj.res$pc.Cov
  L <- obj.res$L


  R <- R1 <- A <- Lambda <- Q <- wt <- vector("list",n.obj)
  for(i in 1:n.obj){
    R[[i]] <- obj.list[[i]]$R
    R1[[i]] <- mat.sqrt(R[[i]])
  }
 

  T <- resampling.test.stat.kins_adj.normal_appx(obj.res,obj.list,Z1,resample)
  Q <- T$Test.Stat

  for(i in 1:n.obj)
	A[[i]] <- obj.list[[i]]$W
  
  resid = NULL

  pv <- matrix(0,ncol = n.obj,nrow = resample)
  for(j in 1:n.obj){
   for(i in 1:resample){
    pv[i,j] <- SKAT:::Get_Davies_PVal(Q[i,j]/2, A[[j]], resid)$p.value
   }
  }

  return(list(p.values = pv,Sigma_Ps = R))
}
