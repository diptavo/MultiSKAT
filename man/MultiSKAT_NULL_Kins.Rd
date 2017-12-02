\name{MultiSKAT_NULL_Kins}
\alias{MultiSKAT_NULL_Kins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct the null MultiSKAT model with kinship adjustment
}
\description{
Compute the null model parameters and residuals for MultiSKAT_Kins. The values and residuals from
here are used to compute the MultiSKAT_Kins test (MultiSKAT with kinship adjustment).
}
\usage{
MultiSKAT_NULL_Kins(y.mat, X, U, D, V_g,V_e)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.mat}{
A numeric matrix of phenotype values with rows as individuals and columns as
phenotypes.
}
  \item{X}{
A numeric matrix of covariate values with rows as individuals and columns as
covariates. The first column should be the intercept (a vector of 1's)
}
  \item{U}{
Spectral vector matrix of the Kinship matrix.
}
  \item{D}{
Spectral values matrix of the Kinship matrix
}
  \item{V_g}{
Co-heritability matrix
}
  \item{V_e}{Residual variance matrix}
}
\details{
Estimation of the coheritability and residual variance matrices should be as accurate as possible, otherwise it can hamper the false positive rates of the corresponding MultiSKAT tests.
}

\value{
This function returns an null object of the class MultiSKAT_NULL_related. This object contains details, values and residuals for computing MultiSKAT test statistic and the null distribution.
}

\author{
Diptavo Dutta
}

\examples{
data(MultiSKAT.Kinship.example.data)
attach(MultiSKAT.Kinship.example.data)

obj.null = MultiSKAT_NULL_Kins(Phenotypes,Cov,U,D,V_g,V_e)


}
\keyword{ Null_model }
