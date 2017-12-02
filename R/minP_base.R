minP_base <- function(obj.res,obj.list,Z1,resample = 200){

  X <- obj.res$cov
  Q <- resampling.pvalues(obj.res,obj.list,Z1,X,resample = 200,Method = "Normal.Appx")$p.values
  n.obj <- length(obj.list);
  n.pheno <- obj.res$n.pheno
  n <- obj.res$n.all

  if(n.obj > 4){
	msg<-sprintf("Combining more than 4 Multi-SKAT tests through minP is not advised!!")
	warning(msg, call. = FALSE)
  }
  pv <- array()
  for(i in 1:n.obj)
	pv[i] <- obj.list[[i]]$p.value

  minP <- min(pv)
  p.minp = n.obj*minP


  ind <-  adj.Lambda(eigen(cov(Q))$values);
  if(length(ind) < n.obj){
	return(list(p.value = p.minp))}else{

  R.mat <- cor(Q,method = "kendall")
  obj.t.Copula = tCopula(R.mat[lower.tri(R.mat)],dim = nrow(R.mat),dispstr = "un");
  obj.norm.Copula = normalCopula(R.mat[lower.tri(R.mat)],dim = nrow(R.mat),dispstr = "un");

  p.minp <- 1-pCopula(1-minP,obj.t.Copula)
  if((p.minp/minP > n.obj)||(p.minp/minP < 1))
	p.minp <-  1-pCopula(1-minP,obj.norm.Copula)
  
  if((p.minp/minP > n.obj)||(p.minp/minP < 1))
	p.minp <-  n.obj*minP

  return(list(p.value = p.minp))}
}

minP_Kins <- function(obj.res,obj.list,Z1,resample = 200){

  X <- obj.res$cov
  Q <- resampling.pvalues.kins_adj(obj.res,obj.list,Z1,resample)$p.values
  n.obj <- length(obj.list);
  n.pheno <- obj.res$n.pheno
  n <- obj.res$n.all
  
  pv <- array()
  for(i in 1:n.obj)
	pv[i] <- obj.list[[i]]$p.value
  

  minP <- min(pv)
  p.minp = n.obj*minP


  ind <-  adj.Lambda(eigen(cov(Q))$values);
  if(length(ind) < n.obj){
	return(list(p.value = p.minp))}else{

  R.mat <- cor(Q,method = "kendall")
 
  obj.t.Copula = tCopula(R.mat[lower.tri(R.mat)],dim = nrow(R.mat),dispstr = "un");
  obj.norm.Copula = normalCopula(R.mat[lower.tri(R.mat)],dim = nrow(R.mat),dispstr = "un");

  p.minp <- 1-pCopula(1-minP,obj.t.Copula)
  if((p.minp/minP > n.obj)||(p.minp/minP < 1))
	p.minp <-  1-pCopula(1-minP,obj.norm.Copula)
  
  if((p.minp/minP > n.obj)||(p.minp/minP < 1))
	p.minp <-  n.obj*minP

  return(list(p.value = p.minp))}
}
