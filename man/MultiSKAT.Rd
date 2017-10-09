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
MultiSKAT(obj.res, Z1, Sigma_p,Is.Common = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.res}{
a NULL object output from MultiSKAT_NULL
}

  \item{Z1}{
a kernelized genotype object output from Genotype.Kernels
}

  \item{Sigma_p}{
a phenotype kernel. Denotes the kernel for the effect sizes of a particular variant in the gene on the battery of phenotypes. 
}
 \item{Is.Common}{
a binary variable indiciating whether a variant has the same effect on all the phenotypes (default=FALSE). 
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

obj.null <- MultiSKAT_NULL_Fast(Phenotypes,Cov)
G.Kernel <- Genotype.Kernels(Genotypes,obj.null)

out <- MultiSKAT(obj.null,G.Kernel,cov(Phenotypes))
out$p.value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
