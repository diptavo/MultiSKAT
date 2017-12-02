MultiSKAT_Kins = function(obj.res,Z,Sigma_p, kernel = "linear.weighted", Is.Common=FALSE, weights.beta=c(1,25), weights = NULL
                            , impute.method = "fixed",r.corr=0,   is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose = TRUE){
  
    if(class(obj.res) != "MultiSKAT_NULL_related"){
	stop("This function should be applied to related samples !!")}else{

  Z1 <- Genotype.Kernels(Z,obj.res, kernel = kernel, Is.Common=Is.Common, weights.beta=weights.beta, weights = weights, impute.method = impute.method,r.corr=r.corr,   is_check_genotype=is_check_genotype, is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,max_maf=max_maf, verbose = verbose)



  R = Sigma_p
  R1 = mat.sqrt(Sigma_p)
  U <- obj.res$U; D <- obj.res$D;V_g <- obj.res$V_g
  
  Z <- Z1
  Z1 <- t(U)%*%Z1
  U1X <- as.matrix(obj.res$U1X)
  
  d <- diag(D)
  m = ncol(Z1)
  n = nrow(Z1)
  n.pheno = obj.res$n.pheno
  q1 = obj.res$q1
  N = obj.res$n.all
  V.item = solve(obj.res$V.item.inv)
  Q = Test.Stat.Mod.Kins_adj(obj.res,Z,Sigma_p,Is.Common = FALSE)$Q
  m1 = Test.Stat.Mod.Kins_adj(obj.res,Z1,Sigma_p,Is.Common = FALSE)$m1
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
    
    ZVinvZ<-ZVinvZ + t(Z2) %*% solve(V.item + d[i]*V_g) %*% Z2
    XVinvZ<-XVinvZ + t(U1X[id,]) %*%solve(V.item + d[i]*V_g) %*% Z2
    
  }
  
  W.1= ZVinvZ - t(XVinvZ) %*% obj.res$XUVU1X1_inv %*% XVinvZ
  # compute p-value using davies method
  out<-SKAT:::Get_Davies_PVal(Q/2, W.1, resid)   
  return(list(Q = Q,p.value = out$p.value,W = W.1,R = Sigma_p))}
}  


Test.Stat.Mod.Kins_adj = function(obj.res,Z1,Sigma_p,Is.Common = FALSE){
  R = Sigma_p
  R1 = mat.sqrt(Sigma_p)
  m = ncol(Z1)
  n = nrow(Z1)
  n.pheno = obj.res$n.pheno
  q1 = obj.res$q1
  N = obj.res$n.all
  V.item.inv = obj.res$V.item.inv
  
  res.mat = matrix(obj.res$kins.adj.res.V, byrow=TRUE, ncol=n.pheno)
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

