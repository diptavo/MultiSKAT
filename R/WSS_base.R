WSS_base <- function(obj.res,obj.list,Z1){
  
  X <- obj.res$cov
  temp1.1<-t(X) %*% Z1
  temp1.2<-solve(t(X)%*%X)
  A1<-t(Z1)%*%Z1 - t(temp1.1) %*% (temp1.2 %*% temp1.1)
  
  V.item.inv <- obj.res$V.item.inv;
  n.obj <- length(obj.list)
  R <- R1 <- A <- Lambda <- Q <- wt <- vector("list",n.obj)
  for(i in 1:n.obj){
    R[[i]] <- obj.list[[i]]$R
    R1[[i]] <- mat.sqrt(R[[i]])
    A[[i]] <- R1[[i]]%*%V.item.inv%*%R1[[i]]
    Lambda[[i]] <- adj.Lambda(as.vector(eigen(A[[i]])$values%*%t(eigen(A1)$values)))
    Q[[i]] <- obj.list[[i]]$Q
    wt[[i]] <- sqrt(1/sum(Lambda[[i]]))
  }
 w <- sum(unlist(wt))
 wt <- unlist(wt)/w
 for(i in 1:n.obj){
   Q[[i]] <- wt[i]*Q[[i]]
   A[[i]] <- wt[i]*A[[i]]
 }
 Q.wss <- sum(unlist(Q))
 A2.wss <- sum(unlist(A)) 
 mat.final.wss <- as.vector(eigen(A2.wss)$values%*%t(eigen(A1)$values))
 
 p1 <- SKAT:::Get_PValue.Lambda(mat.final.wss,Q.wss)$p.value
 return(list(Q = Q.wss,p.value = p1,Lambda = mat.final.wss))
    
}

WSS_Kins <- function(obj.res,obj.list,Z1){

  X <- obj.res$cov	
  V.item.inv <- obj.res$V.item.inv;
  n.obj <- length(obj.list)
  R <- R1 <- W <- Lambda <- Q <- wt <- vector("list",n.obj)
  for(i in 1:n.obj){
    R[[i]] <- obj.list[[i]]$R
    R1[[i]] <- mat.sqrt(R[[i]])
    W[[i]] <- obj.list[[i]]$W
    Q[[i]] <- obj.list[[i]]$Q
    Lambda[[i]] <- sqrt(1/sum(SKAT:::Get_Lambda(W[[i]])))
  }
 w <- sum(unlist(Lambda))
 wt <- unlist(Lambda)/w

 W.mat <- matrix(0,ncol = dim(W[[1]])[1],nrow = dim(W[[1]])[2])
 for(i in 1:n.obj){
   Q[[i]] <- wt[i]*Q[[i]]
   W.mat <- W.mat + wt[i]*W[[i]]
 }
 Q.wss <- sum(unlist(Q))
 A2.wss <- sum(unlist(W)) 
 p <- SKAT:::Get_Davies_PVal(Q.wss/2,W.mat,NULL)$p.value
 
 return(list(Q = Q.wss,p.value = p, W = W.mat,weights = wt))

}


adj.Lambda <- function(lambda1){
    IDX1 <- which(lambda1 >= 0)
    IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
    if (length(IDX2) == 0) {
        return(0)
    }else{
    lambda <- lambda1[IDX2]
    return(lambda)}
}


