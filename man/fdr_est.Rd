\name{fdr_est}
\alias{fdr_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to compute the estimated pointwise FDR for every observed p-max}
\description{
A function to compute the estimated pointwise FDR based on the proposed joint significance mixture null method (JS-mixture).
}
\usage{
fdr_est(alpha00, alpha01, alpha10, alpha1, alpha2, input_pvalues, exact = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha00}{
A numeric number represents the proportion of null \eqn{H_{00}}
}
  \item{alpha01}{
A numeric number represents the proportion of null \eqn{H_{01}}
}
  \item{alpha10}{
A numeric number represents the proportion of null \eqn{H_{10}}
}
  \item{alpha1}{
A numeric number represents the proportion of null alpha=0
}
  \item{alpha2}{
A numeric number represents the proportion of null beta=0
}
  \item{input_pvalues}{
A matrix contains two columns of p-values for candidate mediators.
Column 1 is the p-value of testing if an exposure  is associated with 
the mediator (alpha!=0). 
Column 2 is the p-value of testing if a mediator is associated with 
the outcome adjusted for the exposure (beta!=0) 
}
  \item{exact}{
The option to choose from two methods. exact=0: approximation without 
estimating the CDFs; exact=1: estimate the CDFs nonparametrically}
}
\details{
A function to estimate the pointwise FDR based on the proposed  method 
to estimate the mixture null distribution. See Dai et al (2020) for details of how to compute quantiles using the approximation method (exact=0) or the exact method (exact=1).}
\value{
The estimated pointwise FDR for p-max
}
\references{
James Y. Dai, Janet L. Stanford, Michael LeBlanc. A multiple-testing 
procedure for high-dimensional mediation hypotheses, 
Journal of the American Statistical Association, 2020, DOI: 10.1080/01621459.2020.1765785. 
}
\author{
James Y. Dai and X. Wang
}

\examples{
data(snp_input)
input_pvalues <- snp_input
#To save time for illustration, we use 10 percent of rows
input_pvalues <- input_pvalues[sample(1:nrow(input_pvalues),size=ceiling(nrow(input_pvalues)/10)),]

nullprop <- null_estimation(input_pvalues)
fdr <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
               nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)

}

\keyword{composite null}
\keyword{joint significance}
\keyword{mediation}

