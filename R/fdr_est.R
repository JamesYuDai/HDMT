fdr_est <-function(alpha00,alpha01,alpha10,alpha1,alpha2,input_pvalues,exact=0){

  ## alpha10,alpha01,alpha00 are estimated three types of null proportions
  ## alpha1 is the marginal null proportion for first p-value
  ## alpha2 is the marginal null proportion for second p-value
  ## input pvalues are two columns of p-values
  ## alpha is the level of FWER to be control at
  ## exact=0 corresponding to the approximation used in section 2.2-2.3 in the paper, the default value for exact is 0
  ## exact=1 corresponding to the exact used in section 2.4 in the paper
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
  efdr1 <- rep(0,nmed)

  if (exact==0) {
   for (i in 1:nmed) {
    fdr11 <-  (pmax[i]*alpha01)/mean(pmax<=pmax[i])
    fdr12 <-  (pmax[i]*alpha10)/mean(pmax<=pmax[i])
    fdr2  <-  (pmax[i]*pmax[i]*alpha00)/mean(pmax<=pmax[i])
    efdr1[i] <- fdr11+fdr12+fdr2
   }
  }
  if (exact==1) {
    #library(fdrtool)
    ish11 <- qvalue(input_pvalues[,1])$qvalue<0.25 & qvalue(input_pvalues[,2])$qvalue<0.25


    out1 <- input_pvalues[!ish11,]
    nmed1 <- nrow(out1)

    nmed  <- nrow(input_pvalues)
    cdf12 <- input_pvalues

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

    orderq1 <- pmax
    orderq2 <- pmax

    gcdf1 <- pmax
    gcdf2 <- pmax
    for (i in 1:length(xknots1)) {
      if (i==1) {
        gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]]
      } else {
        if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
          temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]]
          gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
        }
      }
    }

    for (i in 1:length(xknots2)) {
      if (i==1) {
        gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]]
      } else {
        if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
          temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]]
          gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
        }
      }
    }


    gcdf1 <- ifelse(gcdf1>1,1,gcdf1)
    gcdf2 <- ifelse(gcdf2>1,1,gcdf2)

    cdf12[,1] <- gcdf1
    cdf12[,2] <- gcdf2

    for (i in 1:nmed) {
      fdr11 <-  (pmax[i]*cdf12[i,2]*alpha01)/mean(pmax<=pmax[i])
      fdr12 <-  (pmax[i]*cdf12[i,1]*alpha10)/mean(pmax<=pmax[i])
      fdr2  <-  (pmax[i]*pmax[i]*alpha00)/mean(pmax<=pmax[i])
      efdr1[i] <- fdr11+fdr12+fdr2
    }
  }

  efdr1.order <- efdr1[order(pmax,decreasing=T)]

  for (i in 2:nmed)  {
    efdr1.order[i] <- min(efdr1.order[i],efdr1.order[i-1])
  }

  efdr1 <- efdr1.order[rank(-pmax)]
  return(efdr=efdr1)
}

