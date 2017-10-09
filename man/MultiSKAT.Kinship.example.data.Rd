\name{MultiSKAT.Kinship.example.data}
\alias{MultiSKAT.Kinship.example.data}
\docType{data}
\title{
Example data for MultiSKAT with related samples
}
\description{
This is a dataset of 5 phenotypes for 500 individuals genotyped across 20 variants. For this dataset, variants at 2, 5, 11 and 20 are set to be causal. Phenotypes 1, 2 and 3 are associated with the causal SNPs. The effect sizes of individual SNPs are proportional to logarithm of their minor allele frequencies.
}
\usage{data("MultiSKAT.Kinship.example.data")}
\format{
MultiSKAT.example.data contains the following objects:

\itemize{
 \item{Phenotypes}{
: A matrix of 5 phenotypes for 500 individuals with individuals as rows and phenotypes as columns. 
}

\item{Genotypes}{
: A matrix of Genotype values for 500 individuals across 20 variants with individuals as rows and variants as columns
}

\item{Cov}{
: A vector of intercept covariates for 500 individuals.}

\item{Kinship}{: A Kinship matrix for the 500 individuals}
\item{V_g}{: A co-heritability matrix for 5 phenotypes}
\item{V_e}{: A residual variance matrix for the 5 phenotypes}
\item{U}{: Spectral vector matrix of the Kinship matrix}
\item{D}{: Spectral values matrix of the Kinship matrix}

}
}

\examples{
data(MultiSKAT.Kinship.example.data)
attach(MultiSKAT.Kinship.example.data)
Phenotypes[1:4,1:5]
}

\keyword{datasets}
