\name{MultiSKAT.example.data}
\alias{MultiSKAT.example.data}
\docType{data}
\title{
Example data for MultiSKAT
}
\description{
This is a dataset of 5 phenotypes for 5000 individuals genotyped across 56 variants. The correlation between any two variant is 0.7. For this dataset, variants at 16, 39 and 51 are set to be causal. Phenotypes 1, 3 and 4 are associated with the causal SNPs. The effect sizes of individual SNPs are proportional to logarithm of their minor allele frequencies.
}
\usage{data("MultiSKAT.example.data")}
\format{
MultiSKAT.example.data contains the following objects:



\itemize{
 \item{Phenotypes}{
: A matrix of 5 phenotypes for 5000 individuals with individuals as rows and phenotypes as columns. 
}

\item{Genotypes}{
: A matrix of Genotype values for 5000 individuals across 56 variants with individuals as rows and variants as columns. 
}

\item{Cov}{
: A vector of intercept covariates for 5000 individuals.}
}

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(MultiSKAT.example.data)
attach(MultiSKAT.example.data)
Phenotypes[1:4,1:5]
}
\keyword{datasets}
