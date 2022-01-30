adjust_quantile <-
  function(alpha00,alpha01,alpha10,alpha1,alpha2,input_pvalues,exact=0){

    ## This function computes the expected quantiles of the mixture null distribution
    ## input_pvalues is the 2-column matrix storing the two sets of p-values
    ## alpha00, alpha01, alpha10 are the estimated proportions of three nulls,
    ## exact=0 corresponding to the approximation used in section 2.2-2.3 of the paper;
    ## exact=1 corresponding to the exact method used in section 2.4 of the paper

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

    #library(fdrtool)
    nmed <- nrow(input_pvalues)

    ## compute the quantiles using the approximation method


    if (exact==0) {

      c <- (-(1:nmed)/nmed)
      b <- alpha10+alpha01
      a <- alpha00
      pnull <- (-b+sqrt(b^2-4*a*c))/(2*a)

    }

    ##  compute the quantiles using the exact method
    if (exact==1) {

      ish11 <- qvalue(input_pvalues[,1])$qvalue<0.25 & qvalue(input_pvalues[,2])$qvalue<0.25

      cdf12 <- input_pvalues
      orderp1 <- input_pvalues[order(input_pvalues[,1]),1]
      orderp2 <- input_pvalues[order(input_pvalues[,2]),2]


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

      gcdf1 <- orderp1
      gcdf2 <- orderp2

      orderq1 <- orderp1
      orderq2 <- orderp2

      difff <- 1
      ite <- 1

      while(abs(difff)>1e-6 & ite<10) {
        #cat(ite,"..")
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

        cdf12[,1] <- ifelse(gcdf1>1,1,gcdf1)
        cdf12[,2] <- ifelse(gcdf2>1,1,gcdf2)



        c <- (-(1:nmed)/nmed)
        b <- alpha10*cdf12[,1]+alpha01*cdf12[,2]
        a <- alpha00
        pnull <- (-b+sqrt(b^2-4*a*c))/(2*a)

        difff <- max(max(orderq1-pnull),max(orderq2-pnull))

        orderq1 <- pnull
        orderq2 <- pnull
        ite <- ite+1
      }
    }
    return(pnull)
  }
