Test.Stat.Mod <-
function(obj.res,Z1,Sigma_p,Is.Common = F){
  
  R = Sigma_p
  R1 = mat.sqrt(Sigma_p)
  m = ncol(Z1)
  n = nrow(Z1)
  n.pheno = obj.res$n.pheno
  q1 = obj.res$q1
  N = obj.res$n.all
  V.item.inv = obj.res$V.item.inv
  
  res.mat = matrix(obj.res$res.V, byrow=TRUE, ncol=n.pheno)
  Q.temp1<-t(res.mat) %*% Z1;
  
  if(Is.Common){
    Q<-sum(colSums(Q.temp1)^2)/2;
    m1=m;
  }
  if(!Is.Common){
    tr <- matrix(0,ncol = m, nrow = n.pheno);
    #R <- Sigma_P;
    #       R = cov(matrix(obj.res$res.V,ncol = obj.res$n.pheno,byrow=T))
    for(i in 1:m)
      tr[,i] <- R1%*%Q.temp1[,i];
    Q <- sum(tr^2);
    m1 <- m*n.pheno
  }
  resid = NULL;
  return(list(Q = Q,m1 = m1))
}
