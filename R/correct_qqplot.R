correct_qqplot <-
  function(pmax,pnull,opt="all") {
    
    ## pmax is the a vector of max of p-values
    ## pnull is the quantiles corresponding to pmax 
    
    
    ## opt is either "all" - plotting all points or "subset" that select 0.5% data points for a gene-wide CpG analysis
    
    nmed <- length(pmax)
    pindex=1:nmed
    if (opt=="subset")
    {
      ## pick a subset of p-values to plot, avoid overcrowded q-q plots###
      mpoint <- ceiling(quantile(1:length(pmax),0.95)/200)*200
      pindex <- c(1,(1:ceiling(quantile(1:length(pmax),0.95)/200))*200,mpoint:length(pmax))
    }
    
    xmax <- max(c(-log(pnull[order(pnull,decreasing=T)],base=10)[pindex],-log10(1/nmed)))
    ymax <- max(c(-log(pmax[order(pmax,decreasing=T)],base=10)[pindex],-log(pmax[order(pmax,decreasing=T)],base=10)[pindex]))
    if (xmax>0.8*ymax & xmax<1.25*ymax)
    {
      xmax <- max(xmax,ymax)
      ymax <- max(xmax,ymax)
    }
    
    plot((-log(pnull[order(pnull,decreasing=T)],base=10))[pindex],(-log(pmax[order(pmax,decreasing=T)],base=10))[pindex],xlab="log base 10 (expected null p-values)", ylab="log base 10 (observed p-values)",col=3,xlim=c(0,xmax),ylim=c(0,ymax))
    points((-log((nmed:1)/nmed,base=10))[pindex],(-log(pmax[order(pmax,decreasing=T)],base=10))[pindex],pch=2,col=2)
    legend(0.1,max(-log(pmax,base=10)),c("Uniform null","Mixture null"),pch=2:1,col=2:3,bty="n")
    
    abline(0,1)
  }
