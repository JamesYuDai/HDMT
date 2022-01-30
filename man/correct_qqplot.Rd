\name{correct_qqplot}
\alias{correct_qqplot}

\title{
A function to draw the corrected quantile-quantile plot for p-max using the expected quantiles}
\description{
A function to draw the corrected quantile-quantile (Q-Q) plots. The corrected quantiles were computed from 
the mixture null distribution (green dots) and the naive quantiles were computed from the uniform distribution 
(red dots).
}
\usage{
correct_qqplot(pmax, pnull, opt="all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pmax}{
The vector for maximum  p-values }  
\item{pnull}{
The quantiles of pmax based on the estimated mixture null distribution
}
  \item{opt}{
Option to draw the plot. opt="all":use all the data points, 
opt="subset": use a subset of the data points, in case there 
are too many points in a genome-wide setting, to avoid drawing an overcrowded Q-Q plot with a prohibitive image size.
}
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
input_pvalues <- input_pvalues[sample(1:nrow(input_pvalues),
                               size=ceiling(nrow(input_pvalues)/10)),]

pmax <- apply(input_pvalues,1,max)
nullprop <- null_estimation(input_pvalues)
pnull1 <- adjust_quantile(nullprop$alpha10,nullprop$alpha01,nullprop$alpha00,
                          nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)
correct_qqplot(pmax,pnull1)
}