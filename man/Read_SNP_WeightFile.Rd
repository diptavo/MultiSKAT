\name{Read_SNP_WeightFile}
\alias{Read_SNP_WeightFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a file with custom weights}
\description{
Read a file with custom weights
}
\usage{
Read_SNP_WeightFile(FileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FileName}{
input file name of a custom weight.}
}
\details{
The file should be a white-space (space or tab) delimitered file with 2 columns: SNP_ID and weight
value.
Please keep in mind that there should be no header!!
}
\value{
Output object has a hash table of SNP IDs and weights.
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

