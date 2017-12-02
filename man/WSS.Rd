\name{WSS}
\alias{WSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Sum of Squares
}
\description{
Computes the Weighted Sum of squares omnibus test for a list of individual MultiSKAT test
}
\usage{
WSS(obj.res, obj.list, Z, weights.beta=c(1,25),
kernel = "linear.weighted",Is.Common=FALSE,
weights = NULL,r.corr = 0, impute.method = "fixed",
is_check_genotype=TRUE, is_dosage = FALSE, 
missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.res}{
a NULL object from MultiSKAT_NULL or MultiSKAT_NULL_Kins
}
  \item{obj.list}{
A list containing output from MultiSKAT or MultiSKAT_Fast function
}
  \item{Z}{
A kernelized genotype matrix
}
  \item{kernel}{
a type of kernel (default= "linear.weighted"). See detail section.

}

  \item{Is.Common}{
a binary variable indiciating whether a variant has the same effect on all the phenotypes (default=FALSE). 
}

  \item{weights.beta}{
a numeric vector of parameters for the beta weights for the weighted kernels. If
you want to use your own weights, please use the "weights" parameter. It will
be ignored if "weights" parameter is not null.
}


  \item{r.corr}{
the \eqn{\rho} parameter for the compound symmetric correlation structure kernels (default=0).
}


  \item{weights}{a numeric vector of weights for the weighted kernels. 
}

  \item{impute.method}{
a method to impute missing genotypes (default= "fixed"). "bestguess" imputes
missing genotypes as most likely values (0,1,2), "random" imputes missing
genotypes by generating binomial(2,p) random variables (p is the MAF), and
"fixed" imputes missing genotypes by assigning the mean genotype values (2p).
}
  \item{is_check_genotype}{
a logical value indicating whether to check the validity of the genotype matrix
Z (default= TRUE). If Z has non-SNP data, please set it FALSE, otherwise you
will get an error message. If it is FALSE and you use weighted kernels, the
weights should be given through the "weights" parameter.
}

  \item{is_dosage}{
a logical value indicating whether the matrix Z is a dosage matrix. If it is TRUE,
the function will ignore "is_check_genotype".
}

  \item{missing_cutoff}{
a cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing
rates higher than the cutoff will be excluded from the analysis.
}

  \item{estimate_MAF}{
a numeric value indicating how to estimate MAFs for the weight calculation and
the missing genotype imputation. If estimate_MAF=1 (default), it uses all
samples to estimate MAFs. If estimate_MAF=2, only samples with non-missing
phenotypes and covariates are used to estimate MAFs}



\item{max_maf}{
a cutoff of the maximum minor allele frequencies (MAF) (default=1, no cutoff).
Any SNPs with MAF > cutoff will be excluded from the analysis}

\item{verbose}{
a binary indicator to display messages (default=TRUE, displays messages)
}

}
\details{
Weighted sum of squares is implemented with weights inversely proportional to the variance.
}

\value{
The function returns a list containing
  \item{Q}{Value of the test statistic}
  \item{p.value}{P-value for the weighted sum of squares test}
  \item{Lambda}{Mixing parameter for the null distribution of the WSS test statistic}
}

\references{
Ionita-Laza, I.*, Lee, S.*, Makarov, V., Buxbaum, J. Lin, X. (2013). 
Sequence kernel association tests for the combined effect of rare and common variants.  
\emph{American Journal of Human Genetics}, 92, 841-853. 
* contributed equally. 
}

\author{
Diptavo Dutta
}

\examples{
data(MultiSKAT.example.data)
attach(MultiSKAT.example.data)

obj.null <- MultiSKAT_NULL(Phenotypes,Cov)
out1 <- MultiSKAT(obj.null,Genotypes,cov(Phenotypes))
out2 <- MultiSKAT(obj.null,Genotypes,diag(5))
out3 <- MultiSKAT(obj.null,Genotypes,matrix(1,ncol = 5,nrow = 5))

### Combining PhC, Het and Hom
obj.list = list(out1,out2,out3)
obj.wss = WSS(obj.null,obj.list,Genotypes)
str(obj.wss)

### Combining PhC and Het
obj.list = list(out2,out1)
obj.wss = WSS(obj.null,obj.list,Genotypes)
str(obj.wss)

}

\keyword{ Omnibus }
\keyword{ WSS }
