\name{resampling.pvalues}
\alias{resampling.pvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get null resampling p-values
}
\description{
Function to calculate the null p-values for a list of MultiSKAT objects
}
\usage{
resampling.pvalues(obj.res, obj.list, Z1, X, 
	resample = 200, Method = "Normal.Appx")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.res}{
a NULL object output from MultiSKAT_NULL
}
  \item{obj.list}{
a list containing output from MultiSKAT or MultiSKAT_Fast function
}
  \item{Z1}{
A kernelized genotype matrix
}

  \item{X}{a matrix of covariates}
  \item{resample}{
Number of resampling iterations to obtain the null p-values(default = 200)
}
  \item{Method}{
A string denoting the method to use for the resampling p-values. Currently the two methods that are available are "Normal.Appx" and "Permutation".
}
}
\details{
The output from this function is used in minP and minPcom tests. "Permutation" method should be used for unrelated individuals only.
}
\value{
This function returns a list containing

  \item{p.values}{A numeric matrix containing the resampling p-values}
  \item{Sigma_Ps}{A list containing the \eqn{\Sigma p} values corresponding to the objects in \emph{obj.list}}

}

\author{
Diptavo Dutta
}

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
