MultiSKAT <-
function(obj.res,Z1,Sigma_p,Is.Common = F){
  
  R = Sigma_p
  R1 = mat.sqrt(Sigma_p)
  m = ncol(Z1)
  n = nrow(Z1)
  n.pheno = obj.res$n.pheno
  q1 = obj.res$q1
  N = obj.res$n.all
  V.item.inv = obj.res$V.item.inv
  Q = Test.Stat.Mod(obj.res,Z1,Sigma_p)$Q
  m1 = Test.Stat.Mod(obj.res,Z1,Sigma_p)$m1
  ZVinvZ = matrix(rep(0, m1*m1), ncol=m1)
  XVinvZ = matrix(rep(0, m1*q1), ncol=m1)
  resid = NULL
  
  for(i in 1:n){
    id<-((i-1)*n.pheno+1):(i*n.pheno)
    
    if(!Is.Common){
      Z2 = t(Z1[i,]) %x% R1
    } else {
      Z2 = t(Z1[i,]) %x% rep(1, n.pheno)
    }
    
    ZVinvZ<-ZVinvZ + t(Z2) %*% V.item.inv %*% Z2
    XVinvZ<-XVinvZ + t(obj.res$X1[id,]) %*% V.item.inv %*% Z2
    
  }
  
  W.1= ZVinvZ - t(XVinvZ) %*% obj.res$XVX1_inv %*% XVinvZ
  # compute p-value using davies method
  out<-SKAT:::Get_Davies_PVal(Q/2, W.1, resid)    
  return(list(Q = Q,p.value = out$p.value,W = W.1,R = Sigma_p))
}
