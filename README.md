# HDMT
### An R package implementing the multiple-testing procedure for high-dimensional mediation analyses 

<br>

Mediation analysis is of rising interest in clinical trials and epidemiology. The advance of high-throughput technologies has made it possible to interrogate molecular phenotypes such 
as gene expression and DNA methylation in a genome-wide fashion, some of which may act as intermediaries of treatment, external exposures and life-style risk factors in the
etiological pathway to diseases or traits. When testing for mediation in high-dimensional studies, properly controlling the type I error rate remains 
a challenge due to the composite null hypothesis. 

<br>

Among existing methods, the joint significance (JS) test is an intersection-union test using the maximum p-value for testing the two parameters, though a naive significance 
rule based on the uniform null p-value distribution (JS-uniform) may yield an overly conservative type I error rate and therefore low power. This is particularly a concern 
for high-dimensional mediation hypotheses for genome-wide molecular intermediaries such as DNA methylation. 

<br>
In this R package we develop a multiple-testing procedure that accurately controls the family-wise error rate (FWER) and the false discovery rate (FDR) for testing 
high-dimensional mediation composite null hypotheses. The core of our procedure is based on estimating the proportions of three types of component null hypotheses and 
deriving the corresponding mixture distribution (JS-mixture) of null p-values.

<br>

<br>

The R package has the following functions
  - `null_estimation()`:  compute the three null proportions 
  - `adjust_quantile()`: compute the adjusted null quantiles for the input p-values, using the mixture null distribution 
  - `fdr_est()`: compute the false discovery rate of the input p-values
  - `fwer_est()`: compute the family-wise error rate of the input p-values
  - `correct_qqplot()`: draw the q-q plots using the corrected null quantiles. 
  
