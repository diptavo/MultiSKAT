\name{MultiSKAT_NULL}
\alias{MultiSKAT_NULL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct and get parameters and residuals from the null model.
}
\description{
Compute the null model parameters and residuals for MultiSKAT. The values and residuals from here are used to compute the MultiSKAT test.
}
\usage{
MultiSKAT_NULL(y.mat, X,is.fast=TRUE,method = "nlminb")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.mat}{
A numeric matrix of phenotype values with rows as individuals and columns as phenotypes. 
}
  \item{X}{
A numeric matrix of covariate values with rows as individuals and columns as covariates. The first column should be the intercept (a vector of 1's).
}
 \item{is.fast}{
A binary variable indicating whether a fast approximation should be used or not (default = TRUE). }
 \item{method}{
Indicates the optimizer to use for the generalized least squares procedure (default = "nlminb". Other option is "optim"). Works only with is.fast = FALSE}
}
\details{
For most of the scenarios with high sample size (more than 2000), is.fast= TRUE should be used. For lower sample size, is.fast=FALSE should be set. There might be cases where the default optimizer 
"nlminb" fails to converge due to numerical issue. Users can switch to another optimzer "optim" in that case by setting method="optim".

}
\value{
This function returns an null object of the class SKAT_Multiphen_Null_PhenCorr. This object contains details, values and residuals for computing MultiSKAT test statistic and the null distribution.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Diptavo Dutta
}
\examples{
data(MultiSKAT.example.data)
attach(MultiSKAT.example.data)

obj.null <- MultiSKAT_NULL(Phenotypes,Cov)

}

\keyword{ Null_model }
