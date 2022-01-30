\name{adjust_quantile}
\alias{adjust_quantile}

\title{
A function to compute the quantiles of the estimated mixture null distribution for pmax using either the approximation
or exact method
%%  ~~function to do ... ~~
}
\description{
A function to generate the quantiles of the estimated three-component mixture null distribution for pmax 
(the maximum of the two p-values  for testing mediation) using either approximation or exact method
}
\usage{
adjust_quantile(alpha00, alpha01, alpha10, alpha1, alpha2, 
input_pvalues, exact = 0)
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
A numeric number represents the proportion of null alpha==0 
(association between exposure and mediator)}
  \item{alpha2}{
A numeric number represents the proportion of null beta==0 
(association between mediator and outcome adjusted for exposure)}
  \item{input_pvalues}{
A matrix contains two columns of p-values for candidate mediators.
Column 1 is the p-value of testing if a exposure is associated with 
the mediator (alpha!=0). 
Column 2 is the p-value of testing if a mediator is associated with 
the outcome adjusted for exposure(beta!=0)  
}
  \item{exact}{
Use the option to choose from two methods. 
exact=0: the approximation method without estimating the CDFs when deriving the mixture null distribution;
exact=1: the exact method to estimate the CDFs nonparametrically when deriving the mixture null distribution}
}
\details{
The is a function to computed the expected quantiles for the observed 
p-max values based on the estimated mixture null distribution. The methodology detail can be found in Dai et al (2020).
}
\value{
A vector contains the expected quantiles of p-values based on the estimated mixture null distribution. See Dai et al (2020) for details of how to compute quantiles using the approximation method (exact=0) or the exact method (exact=1).
}
\references{
James Y. Dai, Janet L. Stanford, Michael LeBlanc. A multiple-testing 
procedure for high-dimensional mediation hypotheses. 
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

pnull <- adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,
                         nullprop$alpha2,input_pvalues,exact=0)
}

\keyword{composite null}
\keyword{joint significance}
\keyword{mediation}