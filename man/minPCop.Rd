\name{minPCop}
\alias{minPCop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum p-value using copula
}
\description{
Function to compute the p-value for the minimum p-value test statistic using Copula}
\usage{
minPCop(obj.res, obj.list, Z1,X,resample = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.res}{
a NULL object from MultiSKAT_NULL or MultiSKAT_NULL_Kins

}
  \item{obj.list}{
A list containing output from MultiSKAT or MultiSKAT_Fast function
}
  \item{Z1}{
A kernelized genotype matrix
}
 \item{X}{A covariate matrix}
\item{resample}{Number of resampling iterations to compute the correlation matrix between the tests}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a numeric value which is the p-value of the minPCop test.}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
