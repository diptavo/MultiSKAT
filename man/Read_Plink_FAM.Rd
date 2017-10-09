\name{Read_Plink_FAM}
\alias{Read_Plink_FAM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Plink FAM and covariates files
}
\description{
Read Plink FAM and covariates files.
}
\usage{
Read_Plink_FAM(Filename, Is.binary = TRUE, flag1 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Filename}{
input file name of plink FAM file
}
  \item{Is.binary}{
if TRUE, the phenotype is binary. If phenotype is continuous, it should be
FALSE
}
  \item{flag1}{
0 represents the default coding of unaffected/affected (1/2) (default=0), and 1
represents 0/1 coding. flag1=1 is the same as -1 flag in plink. Please see the
plink manual.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A dataframe of Family ID (FID), Individual ID (IID), Paternal ID (PID), Maternal ID(MID), Sex,
and Phenotype.
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
