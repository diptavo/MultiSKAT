\name{MultiSKAT_Fast}
\alias{MultiSKAT_Fast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Faster version of Sequence kernel association test for Multiple Phenotypes
}
\description{
A computationally faster implementation for the MultiSKAT function. Tests for association between a set of SNPs/genes and a battery of continuous correlated phenotypes using kernel regression framework using eigen value calculations.
}
\usage{
MultiSKAT_Fast(obj.res, Z1, Sigma_p, X)
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
  \item{X}{Covariate matrix for the test}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a list containing
  \item{Q }{Value of the test statistic}
  \item{p.value }{P-value of MultiSKAT test}
  \item{Lambda}{Numeric vector of eigen values of the phenotype adjusted regional information. These values are used as the mixing parameters for the null distribution}
  \item{R}{The phenotype kernel that was used for the test}
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
data(MultiSKAT.example.data)
attach(MultiSKAT.example.data)

obj.null <- MultiSKAT_NULL(Phenotypes,Cov)
G.Kernel <- Genotype.Kernels(Genotypes,obj.null)

out <- MultiSKAT_Fast(obj.null,G.Kernel,cov(Phenotypes),Cov)
out$p.value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
