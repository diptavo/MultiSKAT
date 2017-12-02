minP <- function(obj.res,obj.list,Z,resample = 200,weights.beta=c(1,25), kernel = "linear.weighted", Is.Common=FALSE,  weights = NULL,r.corr = 0,
                            impute.method = "fixed",   is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf = 1,verbose = TRUE){

	  Z1 <- Genotype.Kernels(Z,obj.res,r.corr = r.corr, kernel = kernel, Is.Common=Is.Common, weights.beta=weights.beta, weights = weights, impute.method = impute.method,   is_check_genotype=is_check_genotype, is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,max_maf=max_maf,verbose = verbose)


	if(class(obj.res) == "MultiSKAT_NULL_unrelated"){
		re <- minP_base(obj.res,obj.list,Z1,resample = resample)}
	else{
		re <- minP_Kins(obj.res,obj.list,Z1,resample = resample)}
	return(re)
}

