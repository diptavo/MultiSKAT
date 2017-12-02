mat.sqrt <- function(A)
{
  ei<-eigen(A)
  d<-ei$values
  d<-(d+abs(d))/2
  d2<-sqrt(d)
  ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(ans)
}

MAF <- function(G){
  mf <- array()
  for(i in 1:ncol(G)){
    mf[i] <- sqrt(sum(G[,i] == 0)/nrow(G));
    if(mf[i] > 0.5) 
      mf[i] = 1- mf[i];
  }
  return(mf)
}

MAC <- function(G){
  mc <- array()
  for(i in 1:ncol(G)){
    mc[i] <- sum(G[,i] == 1) + 2*min(sum(G[,i] == 2),sum(G[,i] == 0));
  }
  return(mc);
}

