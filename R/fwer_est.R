fwer_est <-function(alpha10,alpha01,alpha00,alpha1,alpha2,input_pvalues,alpha=0.05,exact=0) {

  ## alpha10,alpha01,alpha00 are estimated null proportions
  ## alpha1 is the estimated marginal null proportion for first p-value
  ## alpha2 is the estimated marginal null proportion for second p-value

  ## input_pvalues is data matrix with two columns of p-values
  ## alpha is the level of FWER to be controlled at, default 0.05
  ## exact=0 corresponding to the approximation used in section 2.2-2.3, the default value is 0
  ## exact=1 corresponding to the exact used in section 2.4
  ## check input

  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) !=2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
    stop("input_pvalues doesn't have valid p-values")

  pmax <- apply(input_pvalues,1,max)
  nmed <- length(pmax)

  ## first compute the approximation fwer cut-off, using it as starting value if approximation method 2 is called
  c <- (-alpha)/nmed
  b <- alpha01+alpha10
  a <- alpha00
  if (exact==0) { fwer_alpha<- (-b+sqrt(b^2-4*a*c))/(2*a) }


  if (exact==1) {
    ish11 <- qvalue(input_pvalues[,1])$qvalue<0.25 & qvalue(input_pvalues[,2])$qvalue<0.25

    out1 <- input_pvalues[!ish11,]
    nmed1 <- nrow(out1)

    xx1 <- c(0,out1[order(out1[,1]),1])
    yy1 <- c(0,seq(1,nmed1,by=1)/nmed1)

    gfit1<- gcmlcm(xx1,yy1,type="lcm")
    xknots1 <- gfit1$x.knots[-1]
    Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)

    xx2 <- c(0,out1[order(out1[,2]),2])
    yy2 <- c(0,seq(1,nmed1,by=1)/nmed1)

    gfit2<- gcmlcm(xx2,yy2,type="lcm")
    xknots2 <- gfit2$x.knots[-1]
    Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)

    alpha1 <- (alpha00+alpha01)/(alpha00+alpha01+alpha10)
    alpha2 <- (alpha00+alpha10)/(alpha00+alpha01+alpha10)

    if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
    if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))

    fwer_alpha<- (-b+sqrt(b^2-4*a*c))/(2*a)
    qfwer <- fwer_alpha


    ite <- 1
    difff <- 1
    while (abs(difff)>1e-6 & ite<10){
      cat(ite,"..")
      if (sum(input_pvalues[,1]<fwer_alpha)<60 & sum(input_pvalues[,1]<fwer_alpha)>15) {
        if (alpha1==1) cdf1 <- 0 else {
          cdf1 <- max(0,(mean(input_pvalues[,1]<fwer_alpha) - alpha1*fwer_alpha)/(1-alpha1))
          cdf1 <- min(cdf1,1)
        }
      }  else{
        if (sum(input_pvalues[,1]<fwer_alpha)<=15) cdf1 <- 1
        if (sum(input_pvalues[,1]<fwer_alpha)>=60) {
          if (fwer_alpha<=xknots1[1]) cdf1 <- Fknots1[1] else {
            for (i in 2:length(xknots1)) {
              if (sum(fwer_alpha>xknots1[i-1] & fwer_alpha<=xknots1[i])>0){
                cdf1 <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(fwer_alpha-xknots1[i-1])
              }
            }
            if (fwer_alpha>xknots1[length(xknots1)]) cdf1 <- 1
          }
        }
      }

      if (sum(input_pvalues[,2]<fwer_alpha)<60 & sum(input_pvalues[,2]<fwer_alpha)>15) {
        if (alpha2==1) cdf2 <- 0 else {
          cdf2 <- max(0,(mean(input_pvalues[,2]<fwer_alpha) - alpha2*fwer_alpha)/(1-alpha2))
          cdf2 <- min(cdf2,1)
        }
      }  else{
        if (sum(input_pvalues[,2]<fwer_alpha)<=15) cdf2 <- 1
        if (sum(input_pvalues[,2]<fwer_alpha)>=60) {
          if (fwer_alpha<=xknots2[1]) cdf2 <- Fknots2[1] else {
            for (i in 2:length(xknots2)) {
              if (sum(fwer_alpha>xknots2[i-1] & fwer_alpha<=xknots2[i])>0){
                cdf2 <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(fwer_alpha-xknots2[i-1])
              }
            }
            if (fwer_alpha>xknots2[length(xknots2)]) cdf2 <- 1
          }
        }
      }
      c <- (-alpha/nmed)
      if (cdf1>1) cdf1 <- 1
      if (cdf2>1) cdf2 <- 1

      b <- alpha10*cdf1+alpha01*cdf2
      a <- alpha00
      fwer_alpha <- (-b+sqrt(b^2-4*a*c))/(2*a)
      difff <- max(qfwer-fwer_alpha)
      qfwer <- fwer_alpha
      ite <- ite+1
    }
  }
  return(fwer_alpha)
}

