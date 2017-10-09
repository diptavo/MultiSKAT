MultiSKAT_NULL<-
function(y.mat, X) {
  
  Corr <- cor(y.mat)
  Cov <- cov(y.mat)
  require(nlme)
  # Check dimensions
  X<-cbind(X)  
  n<-nrow(y.mat)
  n.pheno<-ncol(y.mat)
  q<-ncol(X)
  
  if(sum(scale(X[,1], scale = FALSE)^2) != 0){
    stop("The first column of X should be an intercept!")
  }
  
  y<-as.vector(t(y.mat))
  X1<- X %x% diag(n.pheno)
  subject.id<-rep(1:n, each=n.pheno)
  
  if( sum(c(cor(y.mat)[lower.tri(cor(y.mat))]) > 1e-07 ) == 0){
    V.item = cov(y.mat)
    V.item.inv<-solve(V.item)
    gls1 <- lm(y ~ X1-1)
    res <- residuals(gls1)  }else{
      gls1<-gls(y ~ X1-1,  correlation = corSymm(form = ~ 1 | subject.id))
      V.item<-getVarCov(gls1)
      V.item.inv<-solve(V.item)
      res<-residuals(gls1)
    }
  
  
  V.item.inv = solve(Cov)
  
  q1<-ncol(X1)
  XVX1<-matrix(rep(0,q1*q1), ncol=q1)
  res.V<-rep(0,n*n.pheno)
  for(i in 1:n){
    id<-((i-1)*n.pheno+1):(i*n.pheno)
    XVX1<-XVX1 + t(X1[id,]) %*% V.item.inv %*% X1[id,]
    
    res.V[id]<-V.item.inv %*% res[id]
    
  }
  XVX1_inv<-solve(XVX1)
  
  L <- princomp(y.mat)$loadings
  pc.mat <- y.mat%*%L
  V_p = cov(pc.mat)
  
  re<-list(res = res, res.V = res.V, V.item.inv = V.item.inv, X1=X1, XVX1_inv=XVX1_inv, n.pheno=n.pheno, q1=q1, n.all=n, id_include=1:n,y.Corr = Corr,y.Cov = Cov,pc.Cov = V_p,L = L)
  
  class(re)<-"SKAT_Multiphen_Null_PhenCorr"
  return (re)
}