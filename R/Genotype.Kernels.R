Genotype.Kernels <-
function(Z,obj.res, kernel = "linear.weighted", Is.Common=FALSE, weights.beta=c(1,25), weights = NULL
                            , impute.method = "fixed",r.corr=0,   is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose = TRUE){
  
#  kernel = "linear"; Is.Common=FALSE; weights.beta=c(0.5,0.5); weights = NULL  ; impute.method = "fixed";r.corr=0;   is_check_genotype=TRUE; is_dosage = FALSE; missing_cutoff=0.15; estimate_MAF=1; max_maf=1;
  
  m = ncol(Z)
  n = nrow(Z)
 
  if(verbose){ print(paste("The region has ",m," variants",sep = ""))}
  if(m < 3 & verbose){
	msg <-sprintf("Rare variant test with < 3 variants is not advisable")
	warning(msg,call. = FALSE)}
  n.rare <- sum(MAF(Z) < 0.01, na.rm = TRUE)
    if(verbose){print(paste("The region has ",n.rare," rare variants",sep = ""))}
  mc <- sum(MAC(Z), na.rm = TRUE)
  if(mc < 5 & verbose){
	msg <-sprintf("Rare variant test with total MAC < 5 is not advisable")
	warning(msg,call. = FALSE)}

  n.pheno = obj.res$n.pheno
  q1 = obj.res$q1
  V.item.inv = obj.res$V.item.inv
  out.z<-SKAT:::SKAT_MAIN_Check_Z(Z, obj.res$n.all, obj.res$id_include, SetID=NULL, weights, weights.beta, impute.method
                                  , is_check_genotype, is_dosage, missing_cutoff,estimate_MAF=estimate_MAF,max_maf=max_maf)
  if(out.z$return ==1){
    out.z$param$n.marker<-m
    return(out.z)
  }
  
  if (kernel == "linear.weighted") {
    Z1 = t(t(out.z$Z.test) * (out.z$weights))
  } else{
    Z1 = out.z$Z.test
  }
  
  
  if(r.corr == 1){
    Z1<-cbind(rowSums(Z1))
    colnames(Z1) <- "Hom"
  } else if(r.corr > 0){
    
    p.m<-dim(Z1)[2]  
    R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
    L<-chol(R.M,pivot=TRUE)
    Z1<- Z1 %*% t(L) 
  }
  
  return(Z1)
}
