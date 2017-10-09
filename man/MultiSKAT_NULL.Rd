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
MultiSKAT_NULL(y.mat, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.mat}{
A numeric matrix of phenotype values with rows as individuals and columns as phenotypes. 
}
  \item{X}{
A numeric matrix of covariate values with rows as individuals and columns as covariates. The first column should be the intercept (a vector of 1's).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns an null object of the class SKAT_Multiphen_Null_PhenCorr. This object contains details, values and residuals for computing MultiSKAT test statistic and the null distribution.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
