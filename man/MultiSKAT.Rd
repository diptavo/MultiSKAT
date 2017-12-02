\name{MultiSKAT}
\alias{MultiSKAT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequence kernel association test for Multiple Phenotypes
}
\description{
Test for association between a set of SNPs/genes and a battery of continuous correlated phenotypes using kernel regression framework.

}
\usage{
MultiSKAT(obj.res, Z, Sigma_p, kernel = "linear.weighted", Is.Common=FALSE, 
weights.beta=c(1,25), weights = NULL, impute.method = "fixed", r.corr=0, 
is_check_genotype=TRUE,is_dosage = FALSE, missing_cutoff=0.15, 
estimate_MAF=1,max_maf=1,verbose =TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.res}{
a NULL object output from MultiSKAT_NULL
}

  \item{Z}{
a genotype matrix corresponding to the region to be tested
}

  \item{Sigma_p}{
a phenotype kernel. Denotes the kernel for the effect sizes of a particular variant in the gene on the set of phenotypes. 
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

  \item{weights}{a numeric vector of weights for the weighted kernels. 
}

  \item{impute.method}{
a method to impute missing genotypes (default= "fixed"). "bestguess" imputes
missing genotypes as most likely values (0,1,2), "random" imputes missing
genotypes by generating binomial(2,p) random variables (p is the MAF), and
"fixed" imputes missing genotypes by assigning the mean genotype values (2p).
}

  \item{r.corr}{
the \eqn{\rho} parameter for the compound symmetric correlation structure kernels (default=0).
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a list containing
  \item{Q }{Value of the test statistic}
  \item{p.value }{P-value of MultiSKAT test}
  \item{W}{Phenotype adjusted regional information. The eigen values of this matrix are used as the mixing parameters for the null distribution}
  \item{R}{The phenotype kernel that was used for the test}
  \item{Lambda}{If MultiSKAT_Fast function is used, instead of phenotype adjusted regional information it returns the eigen values of the same, which are used for p-value calculation.}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(MultiSKAT.example.data)
attach(MultiSKAT.example.data)

obj.null <- MultiSKAT_NULL(Phenotypes,Cov)

out <- MultiSKAT(obj.null,Genotypes,cov(Phenotypes))
out$p.value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MultiSKAT }
