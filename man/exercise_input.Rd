\name{exercise_input}
\alias{exercise_input}
\docType{data}
\title{
An example dataset to demonstrate the usage of 'HDMT'
}
\description{
This example dataset was included to assess the mediation role of DNA methylation in the effect of exercise 
on prostate cancer progression in a 
Seattle-based cohort of patients diagnosed with clinically localized 
PCa. The entire data set contains two sets of p-values from genome-wide testing of 450K CpG sites. Due to space limit, 
a subset (10 percent) of the full dataset  is included in the package for illustration.

The dataset is a matrix containing two columns of p-values for 
candidate mediators. Column 1 contains the p-values for testing if an exposure is associated with 
the mediator (alpha!=0). Column 2 contains the p-value for testing if a mediator is associated with 
the outcome after adjusted for the exposure (beta!=0).
}
\usage{data("exercise_input")}
\format{
  The format of exercise_input is:
 num [1:47900, 1:2] 0.4966344 0.1048730 0.1005355 0.4946623 ...
}

\references{

James Y. Dai, Janet L. Stanford, Michael LeBlanc. A multiple-testing 
procedure for high-dimensional mediation hypotheses, 
Journal of the American Statistical Association, 2020, DOI: 10.1080/01621459.2020.1765785.  

}
\examples{
data(exercise_input)
dim(exercise_input)

}
\keyword{dataset}
