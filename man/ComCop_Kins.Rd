\name{ComCop_Kins}
\alias{ComCop_Kins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the p-value for ComCop test
}
\description{
This function combines the choices for \eqn{\Sigma G} and \eqn{\Sigma p} kernel in a single grand omnibus test with Kinship correction. The p-values are obtained by the minimum p-value statistic using copula.
}
\usage{ComCop_Kins(obj.res,Sigma_Ps,Z,X,resample = 200,
weights.beta=c(1,25),kernel = "linear", Is.Common=FALSE, 
weights = NULL, impute.method = "fixed",
is_check_genotype=TRUE, is_dosage = FALSE, 
missing_cutoff=0.15, estimate_MAF=1)
}
\arguments{
  \item{obj.res}{
a NULL object output from MultiSKAT_NULL or MultiSKAT_NULL_Kins
}
  \item{Sigma_Ps}{
A list of \eqn{\Sigma p} kernels to be combined
}
  \item{Z}{
Genotype matrix
}
  \item{X}{
Covariate matrix
}
  \item{resample}{
Number of resampling iterations to obtain the null p-values(default = 200)
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

}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
}
\keyword{datasets}
