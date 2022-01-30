\name{snp_input}
\alias{snp_input}
\docType{data}
\title{
An example dataset to demonstrate the usage of 'HDMT'
}
\description{
This example dataset is included in 'HDMT' to assess the mediation role of DNA methylation 
in genetic regulation of gene expression in primary prostate 
cancer (PCa) samples from The Cancer Genome Atlas (TCGA) with risk SNPs 
as the exposure.

The dataset is a matrix containing two columns of p-values for 
candidate mediators. Column 1 contains the p-values for testing if an exposure is associated with 
the mediator (alpha!=0). Column 2 contains the p-value for testing if a mediator is associated with 
the outcome after adjusted for the exposure (beta!=0).

}
\usage{data("snp_input")}
\format{
  The format of snp_input is:
 num [1:69602, 1:2] 0.106 0.999 0.101 0.173 0.89 ...
}

\references{
James Y. Dai, Janet L. Stanford, Michael LeBlanc. A multiple-testing 
procedure for high-dimensional mediation hypotheses, 
Journal of the American Statistical Association, 2020, DOI: 10.1080/01621459.2020.1765785. 
}
\examples{
data(snp_input)
dim(snp_input)

}
\keyword{dataset}
