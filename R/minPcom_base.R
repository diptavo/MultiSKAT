minPcom_base <- function(obj.res,Sigma_Ps,Z,resample = 200,weights.beta=c(1,25), kernel = "linear.weighted", Is.Common=FALSE,  weights = NULL
                            , impute.method = "fixed", is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose = TRUE){


	Z1 <- Genotype.Kernels(Z,obj.res,r.corr = 0,kernel,weights.beta, Is.Common=FALSE,  weights, impute.method , is_check_genotype, is_dosage, missing_cutoff, estimate_MAF,max_maf,verbose = verbose)   
	Z2 <- Genotype.Kernels(Z,obj.res,r.corr = 1,kernel,weights.beta, Is.Common=FALSE,  weights, impute.method , is_check_genotype, is_dosage, missing_cutoff, estimate_MAF,max_maf,verbose = verbose)   

	X <- obj.res$cov;	
	n.obj <- length(Sigma_Ps);
  	n.pheno <- obj.res$n.pheno
  	n <- obj.res$n.all

	pv <- matrix(0,ncol = 2,nrow = n.obj)
	obj.list1 <- obj.list2 <- vector("list",n.obj)
	for(i in 1:n.obj){
		obj.list1[[i]] <- MultiSKAT(obj.res,Z,Sigma_Ps[[i]],Is.Common = F,r.corr = 0,verbose = verbose); pv[i,1] <- obj.list1[[i]]$p.value
		obj.list2[[i]] <- MultiSKAT(obj.res,Z,Sigma_Ps[[i]],Is.Common = F,r.corr = 1,verbose = verbose); pv[i,2] <- obj.list2[[i]]$p.value
	}
	
	p.res.SKAT <- resampling.pvalues(obj.res,obj.list1,Z1,X,resample,Method = "Normal.Appx")$p.values
   	p.res.Burden <- resampling.pvalues(obj.res,obj.list2,Z2,X,resample,Method = "Normal.Appx")$p.values

	
	R.mat <- cor(cbind(p.res.SKAT,p.res.Burden),method = "kendall")
	minP <- min(pv)
	obj.t.Copula = tCopula(R.mat[lower.tri(R.mat)],dim = nrow(R.mat),dispstr = "un");
  	obj.norm.Copula = normalCopula(R.mat[lower.tri(R.mat)],dim = nrow(R.mat),dispstr = "un");
	
	p.minp <- 1-pCopula(1-minP,obj.t.Copula)
	if((p.minp/minP > 2*n.obj)||(p.minp/minP < 1))
		p.minp <-  1-pCopula(1-minP,obj.norm.Copula)
  
	if((p.minp/minP > 2*n.obj)||(p.minp/minP < 1))
		p.minp <-  2*n.obj*minP
	
	return(list(p.value = p.minp,obj.SKAT = obj.list1,obj.Burden = obj.list2))
}

