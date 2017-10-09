MultiSKAT_Fast <-
function(obj.res,Z1,Sigma_p, X){
  R = Sigma_p
  R1 = mat.sqrt(Sigma_p)
  m = ncol(Z1)
  n = nrow(Z1)
  n.pheno = obj.res$n.pheno
  q1 = obj.res$q1
  N = obj.res$n.all
  
  V.item.inv = obj.res$V.item.inv
  Q = Test.Stat.Mod(obj.res,Z1,Sigma_p)$Q
  
  
  temp1.1<-t(X) %*% Z1
  temp1.2<-solve(t(X)%*%X)
  
  A1<-t(Z1)%*%Z1 - t(temp1.1) %*% (temp1.2 %*% temp1.1)
  #A1 = t(Z1)%*%obj.res$P%*%Z1
  A2 = R1%*%V.item.inv%*%R1
  
  Lambda <- as.vector(eigen(A2)$values%*%t(eigen(A1)$values))
  p1 <- SKAT:::Get_PValue.Lambda(Lambda,Q)$p.value
  return(list(Q = Q,p.value = p1,Lambda = Lambda,R = Sigma_p))  
  
}
