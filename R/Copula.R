resampling.test.stat.normal_appx <- function(obj.res,obj.list,Z1,resample = 200){
  
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
  
  res.mat = matrix(obj.res$res, byrow=TRUE, ncol=n.pheno)
  res <- obj.res$res
  V.item.inv = obj.res$V.item.inv
  V.half <- mat.sqrt(V.item.inv)
  
  res.mat.1 = matrix(obj.res$res.V, byrow=TRUE, ncol=n.pheno)
  Q1 <- Q2 <- Q3 <- Q4 <- array()
  Qs <- matrix(0,ncol = n.obj,nrow = resample)

  for(j in 1:resample){
    
    
    res.1 <- scale(as.matrix(rnorm(n*n.pheno),ncol = n.pheno),center = T,scale = F)
    res.1 <- as.vector(t(res.1))
    
    res.Vhalf<-rep(0,n*n.pheno)
    for(id1 in 1:n){
      idx<-((id1-1)*n.pheno+1):(id1*n.pheno)
      
      res.Vhalf[idx]<-V.half %*% res.1[idx]
      
    }
    
    res.mat.2 = scale(matrix(res.Vhalf, byrow=TRUE, ncol=n.pheno),center = T,scale = F)
    
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


resampling.pvalues <- function(obj.res,obj.list,Z1,X,resample = 200,Method = "Normal.Appx"){
  
  if(Method == "Permutation")
	T <- resampling.test.stat.permutation(obj.res,obj.list,Z1,resample)

  T <- resampling.test.stat.normal_appx(obj.res,obj.list,Z1,resample)
  Q <- T$Test.Stat

  V_y <- obj.res$y.Cov
  V_p <- obj.res$pc.Cov
  L <- obj.res$L
  
  n.obj <- length(obj.list)
  
  n.pheno <- obj.res$n.pheno
  n <- obj.res$n.all
  V.item.inv = obj.res$V.item.inv
  m <- ncol(Z1)
  
  R <- R1 <- A <- Lambda <-  wt <- vector("list",n.obj)
  for(i in 1:n.obj){
    R[[i]] <- obj.list[[i]]$R
    R1[[i]] <- mat.sqrt(R[[i]])
  }

  pv <- matrix(0,ncol = n.obj,nrow = resample) 
  for(i in 1:n.obj)
    pv[[i]] <- array()
  
  temp1.1<-t(X) %*% Z1
  temp1.2<-solve(t(X)%*%X)
  
  A1<-t(Z1)%*%Z1 - t(temp1.1) %*% (temp1.2 %*% temp1.1)

  for(i in 1:n.obj){
    A[[i]] <- R1[[i]]%*%V.item.inv%*%R1[[i]]
    Lambda[[i]] <- as.vector(eigen(A[[i]])$values%*%t(eigen(A1)$values))
  }  
  
  for(j in 1:n.obj){
	for(i in 1:resample)
	    pv[i,j] <- SKAT:::Get_PValue.Lambda(Lambda[[j]],Q[i,j])$p.value
  }   
  
  return(list(p.values = pv,Sigma_Ps = R))
}

resampling.test.stat.permutation <- function(obj.res,obj.list,Z1,resample = 200){
  
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
  
  res.mat = matrix(obj.res$res, byrow=TRUE, ncol=n.pheno)
  res <- obj.res$res
  V.item.inv = obj.res$V.item.inv
  V.half <- mat.sqrt(V.item.inv)
  
  res.mat.1 = matrix(obj.res$res.V, byrow=TRUE, ncol=n.pheno)
  Q1 <- Q2 <- Q3 <- Q4 <- array()
  Qs <- matrix(0,ncol = n.obj,nrow = resample)

  Q1 <- Q2 <- Q3 <- Q4 <- array()
  
  for(j in 1:resample){
    
    
    s <- sample(1:n,n,replace = F)
    res.mat.2 <- res.mat.1[s,]
    
    
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

