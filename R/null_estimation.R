

null_estimation <- function(input_pvalues)
{
  ## updated function that automatically choose best lambda that result in better behave q-q plot

  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) != 2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues)) < nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues),
                                 ]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues) <
      1)
    stop("input_pvalues doesn't have valid p-values")

  ### first identify features that may come from H11 (alternatives for both hypotheses) ###
  #library(qvalue)
  #ish11 <- qvalue(input_pvalues[,1])$qvalue<0.25 & qvalue(input_pvalues[,2])$qvalue<0.25

  pcut <- seq(0.1, 0.8, 0.1)
  frac1 <- rep(0, 8)
  frac2 <- rep(0, 8)
  frac12 <- rep(0, 8)
  for (i in 1:8) {
    frac1[i] <- mean(input_pvalues[, 1] >= pcut[i])/(1 - pcut[i])
    frac2[i] <- mean(input_pvalues[, 2] >= pcut[i])/(1 - pcut[i])
    frac12[i] <- mean(input_pvalues[, 2] >= pcut[i] & input_pvalues[,1] >= pcut[i])/(1 - pcut[i])^2
  }

  alphaout <- matrix(0,4,5)
  ll <- 1
  qqslope <- rep(0,4)
  for (lambda in c(0.5,0.6,0.7,0.8)) {
    alpha00 <- min(frac12[pcut >= lambda][1], 1)
    if (ks.test(input_pvalues[, 1], "punif", 0, 1, alternative = "greater")$p > 0.05)
      alpha1 <- 1 else alpha1 <- min(frac1[pcut >= lambda][1], 1)
      if (ks.test(input_pvalues[, 2], "punif", 0, 1, alternative = "greater")$p > 0.05)
        alpha2 <- 1 else alpha2 <- min(frac2[pcut >= lambda][1], 1)
        if (alpha00 == 1) {
          alpha01 <- 0
          alpha10 <- 0
          alpha11 <- 0
        } else {
          if (alpha1 == 1 & alpha2 == 1) {
            alpha01 <- 0
            alpha10 <- 0
            alpha11 <- 0
            alpha00 <- 1
          }
          if (alpha1 == 1 & alpha2 != 1) {
            alpha10 <- 0
            alpha11 <- 0
            alpha01 <- alpha1 - alpha00
            alpha01 <- max(0, alpha01)
            alpha00 <- 1 - alpha01
          }
          if (alpha1 != 1 & alpha2 == 1) {
            alpha01 <- 0
            alpha11 <- 0
            alpha10 <- alpha2 - alpha00
            alpha10 <- max(0, alpha10)
            alpha00 <- 1 - alpha10
          }
          if (alpha1 != 1 & alpha2 != 1) {
            alpha10 <- alpha2 - alpha00
            alpha10 <- max(0, alpha10)
            alpha01 <- alpha1 - alpha00
            alpha01 <- max(0, alpha01)
            if ((1 - alpha00 - alpha01 - alpha10) < 0) {
              alpha11 <- 0
              alpha10 <- 1 - alpha1
              alpha01 <- 1 - alpha2
              alpha00 <- 1 - alpha10 - alpha01
            }
            else {
              alpha11 <- 1 - alpha00 - alpha01 - alpha10
            }
          }
        }

        pmax <- apply(input_pvalues,1,max)
        pmax <- pmax[order(pmax)]
        nnulls <- sum(pmax>0.8)
        nmed <- nrow(input_pvalues)
        pexp <- rep(0,nnulls)
        for (i in 1:nmed) {
          c <- (-i/nmed)
          b <- alpha01+alpha10
          a <- 1-b
          pexp[i] <- (-b+sqrt(b^2-4*a*c))/(2*a)
        }
        xx <- -log(pexp[(nmed-nnulls+1):nmed],base=10)
        yy <- -log(pmax[(nmed-nnulls+1):nmed],base=10)
        fit1 <- lm(yy~xx-1)

        qqslope[ll]<- fit1$coef[1]
        alphaout[ll,1] <- alpha10
        alphaout[ll,2] <- alpha01
        alphaout[ll,3] <- alpha00
        alphaout[ll,4] <- alpha1
        alphaout[ll,5] <- alpha2

        ll <- ll+1

  }

  bestslope <- which.min(qqslope)
  alpha.null <- list(alpha10 = alphaout[bestslope,1], alpha01 = alphaout[bestslope,2],alpha00 = alphaout[bestslope,3], alpha1 = alphaout[bestslope,4], alpha2 = alphaout[bestslope,5])
  return(alpha.null)
}

