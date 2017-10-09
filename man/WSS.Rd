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
WSS(obj.res, obj.list, Z1, X)
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
  \item{X}{
A matrix of covariates
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a list containing
  \item{Q}{Value of the test statistic}
  \item{p.value}{P-value for the weighted sum of squares test}
  \item{Lambda}{Mixing parameter for the null distribution of the WSS test statistic}
%% ...
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
