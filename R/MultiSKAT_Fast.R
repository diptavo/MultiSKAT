MultiSKAT_Fast <-
function(obj.res,Z,Sigma_p, kernel = "linear.weighted", Is.Common=FALSE, weights.beta=c(1,25), weights = NULL
                            , impute.method = "fixed",r.corr=0,   is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose = TRUE){

  if(class(obj.res) != "MultiSKAT_NULL_unrelated"){
	stop("This function is for independent samples !!")}else{

  Z1 <- Genotype.Kernels(Z,obj.res, kernel = kernel, Is.Common=Is.Common, weights.beta=weights.beta, weights = weights, impute.method = impute.method,r.corr=r.corr,   is_check_genotype=is_check_genotype, is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,max_maf=max_maf,verbose = verbose)


  X <- obj.res$cov
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
  return(list(Q = Q,p.value = p1,Lambda = Lambda,R = Sigma_p))  }
  
}
