\name{Genotype.Kernels}
\alias{Genotype.Kernels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernelize the Genotype matrix for the region to be tested.
}
\description{
This function kernelizes the genotype matrix
}
\usage{
Genotype.Kernels(Z, obj.res, kernel = "linear.weighted", Is.Common = FALSE,
weights.beta = c(1,25), weights = NULL, impute.method = "fixed", 
r.corr = 0, is_check_genotype = TRUE, is_dosage = FALSE, 
missing_cutoff = 0.15, estimate_MAF = 1,max_maf=1,verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
a numeric genotype matrix with each row as a different individual and each
column as a separate gene/snp. Each genotype should be coded as 0, 1, 2, and
9 (or NA) for AA, Aa, aa, and missing, where A is a major allele and a is a
minor allele. Missing genotypes will be imputed by the simple Hardy-Weinberg
equilibrium (HWE) based imputation.}

  \item{obj.res}{
an output object of the MultiSKAT_NULL function.}

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
There are 6 types of pre-specified kernels:
"linear", "linear.weighted", "IBS", "IBS.weighted", "quadratic" and "2wayIX".
Among them, "2wayIX" is a product kernel consisting of main effects and SNP-SNP interaction terms. 
If users want to use dosage values instead of genotypes, set is_dosage=TRUE. 

The r.corr represents a \eqn{\rho} parameter of the unified test, 
\eqn{Q_{\rho} = (1-\rho) Q_S + \rho Q_B}, where \eqn{Q_S} is a SKAT test statistic, 
and \eqn{Q_B} is a weighted burden test statistic. 
Therefore, \eqn{\rho=0} results in the original weighted linear kernel SKAT, 
and \eqn{\rho=1} results in the weighted burden test (default: \eqn{\rho=0}).

If users want to silent the messages from the function, set verbose=FALSE. 

By default, SKAT uses impute.method="fixed" that imputes missing genotypes as the mean genotype values (2p). 
When variates are very rare and missing rates between cases and controls are highly unbalanced, 
impute.method="fixed" can yield inflated type I error rate. In this case, we recommend to use
impute.method="bestguess", which does not suffer the same problem.

}
\value{
This function returns the Kernelized Genotype matrix.
}

\author{
Seunuggeun Lee
}

\references{
Lee, S., Wu, M. C., and Lin, X. (2012)  Optimal tests for rare variant effects in sequencing association studies. \emph{Biostatistics}, 13, 762-775.

Wu, M. C.*, Lee, S.*, Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011)  Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). \emph{American Journal of Human Genetics}, 89, 82-93. \\
* contributed equally. 

}

\examples{
data(MultiSKAT.example.data)
attach(MultiSKAT.example.data)

obj.null <- MultiSKAT_NULL(Phenotypes,Cov)
G.Kernel <- Genotype.Kernels(Genotypes,obj.null)


}
\keyword{ Genotype }
