MultiSKAT <- function(obj.res,Z,Sigma_p, kernel = "linear.weighted", Is.Common=FALSE, weights.beta=c(1,25), weights = NULL
                            , impute.method = "fixed",r.corr=0,   is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose =TRUE){

 if(class(obj.res) == "MultiSKAT_NULL_unrelated"){
	re <- MultiSKAT_Fast(obj.res, Z,Sigma_p = Sigma_p, kernel = kernel, Is.Common=Is.Common, weights.beta=weights.beta, weights = weights, impute.method = impute.method,r.corr=r.corr, is_check_genotype=is_check_genotype, is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,max_maf=max_maf,verbose = verbose)}else{
	re <- MultiSKAT_Kins(obj.res,Z, Sigma_p = Sigma_p, kernel = kernel, Is.Common=Is.Common, weights.beta=weights.beta, weights = weights, impute.method = impute.method,r.corr=r.corr,   is_check_genotype=is_check_genotype, is_dosage = is_dosage, missing_cutoff=missing_cutoff, estimate_MAF=estimate_MAF,max_maf=max_maf,verbose = verbose)}
	class(re) <- "MultiSKAT_Object"
	return(re)
}
