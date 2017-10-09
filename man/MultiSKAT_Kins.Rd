\name{MultiSKAT_Kins}
\alias{MultiSKAT_Kins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequence kernel association test for Multiple Phenotypes with Kinship adjustment
}
\description{
Test for association with kinship adjustment between a set of SNPs/genes and a battery of continuous correlated phenotypes using kernel regression framework.
}
\usage{
MultiSKAT_Kins(obj.res, Z1, Sigma_p, Is.Common = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.res}{
a NULL object output from SKAT_Null_Multiphen_Phencorr_Kins

}
  \item{Z1}{
a kernelized genotype object output from Genotype.Kernels
}
  \item{Sigma_p}{
a phenotype kernel. Denotes the kernel for the effect sizes of a particular variant
in the gene on the battery of phenotypes.
}
\item{Is.Common}{a binary variable indiciating whether a variant has the same effect on all the
phenotypes (default=FALSE).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
