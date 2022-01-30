\name{null_estimation}
\alias{null_estimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to estimate the proportions of the three component nulls
}
\description{
This is a function to estimate the proportions of the three component
nulls involved in mediation testing. We developed a three component-mixture model method to estimate the proportions of nulls and provide much more accurate control of the family-wise error rate (FWER) and the false discovery rate (FDR), when  compared to the standard approach using the uniform null distribution.
}
\usage{
null_estimation(input_pvalues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_pvalues}{
A matrix contains two columns of p-values for candidate mediators.
Column 1 is the p-value of testing if the exposure is associated with 
the candidate mediator (alpha!=0). 
Column 2 is the p-value of testing if the candidate mediator is associated with 
the outcome adjusted for the exposure (beta!=0).
}
}
\details{
A function to estimate the proportions of the three types 
of component null hypotheses:

\eqn{H_{00}}: alpha==0 and beta==0

\eqn{H_{01}}: alpha==0 and beta!=0

\eqn{H_{10}}: alpha!=0 and beta==0

The methodology detail can be found in Dai et al (2020).
}
\value{
A list contains five elements.
\item{alpha00}{A numeric number represents the proportion of null \eqn{H_{00}}}
\item{alpha01}{A numeric number represents the proportion of null \eqn{H_{01}}}
\item{alpha10}{A numeric number represents the proportion of null \eqn{H_{10}}}
\item{alpha1}{A numeric number represents the proportion of null alpha==0}
\item{alpha2}{A numeric number represents the proportion of null beta==0}
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
#To save computing time for illustration, we use 10 percent of rows (p-values)
input_pvalues <- input_pvalues[sample(1:nrow(input_pvalues),
size <- ceiling(nrow(input_pvalues)/10)),]

nullprop <- null_estimation(input_pvalues)
}

\keyword{composite null}
\keyword{joint significance}
\keyword{mediation}
