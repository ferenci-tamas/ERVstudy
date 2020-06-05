prepanel.Dotplot <- function(x, y, horizontal, ...)
{
  if ( horizontal ) {
    xlim <- range(x, attr(x,'other'), na.rm=TRUE)
    ylim <- range(as.numeric(y), na.rm=TRUE)  ## as.numeric 25nov02
  } else {
    xlim <- range(as.numeric(x), na.rm=TRUE)
    ylim <- range(y, attr(y,'other'), na.rm=TRUE)
  }
  list(xlim=xlim, ylim=ylim) #, dx=diff(xlim), dy=diff(ylim))
}

panel.Dotplot <- function(x, y, groups = NULL,
                          pch  = dot.symbol$pch, 
                          col  = dot.symbol$col, cex = dot.symbol$cex, 
                          font = dot.symbol$font, abline, horizontal, ...)
{
  gfun <- ordGridFun(TRUE) ## see Misc.s
  segmnts <- gfun$segments
  y <- as.numeric(y)
  
  gp <- length(groups)
  dot.symbol <- trellis.par.get(if(gp)'superpose.symbol'
                                else 'dot.symbol')
  
  dot.line   <- trellis.par.get('dot.line')
  plot.line  <- trellis.par.get(if(gp)'superpose.line'
                                else 'plot.line')
  
  if ( horizontal )
    gfun$abline(h = unique(y), lwd=dot.line$lwd, lty=dot.line$lty, 
                col=dot.line$col)
  else
    gfun$abline(v = unique(x), lwd=dot.line$lwd, lty=dot.line$lty, 
                col=dot.line$col)
  if(!missing(abline))
  {
    if(length(names(abline))) do.call("panel.abline", abline)
    else for(i in 1:length(abline)) do.call("panel.abline", abline[[i]])
  }
  
  if ( horizontal ) {
    other <- attr(x,'other')
    x <- unclass(x)
    attr(x,'other') <- NULL
  } else {
    other <- attr(y,'other')
    y <- unclass(y)
    attr(y,'other') <- NULL
  }
  
  if(length(other))
  {
    nc <- ncol(other)
    if ( horizontal ) {
      segmnts(other[,1], y, other[,nc], y, lwd=plot.line$lwd[1],
              lty=plot.line$lty[1], col=plot.line$col[1])
      if(nc==4)
      {
        segmnts(other[,2], y, other[,3], y, lwd=2*plot.line$lwd[1],
                lty=plot.line$lty[1], col=plot.line$col[1])
        gfun$points(other[,2], y, pch=3, cex=cex, col=col, font=font)
        gfun$points(other[,3], y, pch=3, cex=cex, col=col, font=font)
      }
    } else {
      segmnts(x, other[,1], x, other[,nc], lwd=plot.line$lwd[1],
              lty=plot.line$lty[1], col=plot.line$col[1])
      if(nc==4)
      {
        segmnts(x, other[,2], x, other[,3], lwd=2*plot.line$lwd[1],
                lty=plot.line$lty[1], col=plot.line$col[1])
        gfun$points(x, other[,2], pch=3, cex=cex, col=col, font=font)
        gfun$points(x, other[,3], pch=3, cex=cex, col=col, font=font)
      }
    }
    
    if(gp) panel.superpose(x, y, groups=as.numeric(groups), pch=pch,
                           col=col, cex=cex, font=font, ...)
    else
      gfun$points(x, y, pch=pch[1], cex=cex, col=col, font=font)
  }
  else
  {
    if(gp) 
      panel.superpose(x, y, groups=as.numeric(groups),
                      pch=pch, col=col, cex=cex,
                      font=font, ...)
    else
      panel.dotplot(x, y, pch=pch, col=col, cex=cex, font=font, horizontal = horizontal, ...)
  }
  if(gp)
  {
    Key <- function(x=0, y=1, lev, cex, col, font, pch, other)
    {
      if(!length(x)) x <- 0.05
      if(!length(y)) y <- 0.95  ## because of formals()
      rlegendg(x, y, legend=lev, cex=cex, col=col, pch=pch, other=other)
      invisible()
    }
    
    lev <- levels(as.factor(groups))
    ng <- length(lev)
    formals(Key) <- list(x=NULL,y=NULL,lev=lev,
                         cex=cex[1:ng], col=col[1:ng],
                         font=font[1:ng], pch=pch[1:ng], other=NULL)
    Hmisc:::.setKey(Key)
  }
}


Dotplot <-
  function (formula, data=sys.frame(sys.parent()),
            groups, subset,
            xlab=NULL, ylab=NULL, ylim=NULL,
            panel=panel.Dotplot, prepanel=prepanel.Dotplot,
            scales=NULL, xscale=NULL, yscale = NULL, ...)
  {
    
    horizontal <- is.factor( data[ , all.vars( formula )[ 1 ] ] )
    
    if ( horizontal ) {
      yvname <- as.character(formula[2])  # tried deparse
      yv <- eval(parse(text=yvname), data)
      if(!length(ylab))
        ylab <- label(yv, units=TRUE, plot=TRUE,
                      default=yvname, grid=TRUE)
      
      if(!length(ylim))
      {
        yother <- attr(yv,'other')
        if(length(yother)) ylim <- range(yv, yother, na.rm=TRUE)
      }
      
      if(is.character(yv)) yv <- factor(yv)
      if(!length(scales) && is.factor(yv))
        scales <- list(y=list(at=1:length(levels(yv)),labels=levels(yv)))
      if(length(xscale)) scales$x <- xscale
      
      xvname <- formula[[3]]
      if(length(xvname)>1 && as.character(xvname[[1]])=='|') 
        xvname <- xvname[[2]]  # ignore conditioning var
      xv <- eval(xvname, data)
      if(!length(xlab)) xlab <- label(xv, units=TRUE, plot=TRUE,
                                      default=as.character(xvname)[1], grid=TRUE)
    } else {
      yvname <- as.character(formula[2])  # tried deparse
      if(length(yvname)>1 && as.character(yvname[[1]])=='|') 
        yvname <- yvname[[2]]  # ignore conditioning var
      yv <- eval(yvname, data)
      if(!length(ylab)) ylab <- label(yv, units=TRUE, plot=TRUE,
                                      default=as.character(yvname)[1], grid=TRUE)
      
      xvname <- formula[[3]]
      xv <- eval(parse(text=xvname), data)
      if(!length(xlab))
        xlab <- label(xv, units=TRUE, plot=TRUE,
                      default=as.character(xvname)[1], grid=TRUE)
      
      if(!length(xlim))
      {
        xother <- attr(xv,'other')
        if(length(xother)) xlim <- range(xv, xother, na.rm=TRUE)
      }
      
      if(is.character(xv)) xv <- factor(xv)
      if(!length(scales) && is.factor(xv))
        scales <- list(x=list(at=1:length(levels(xv)),labels=levels(xv)))
      if(length(yscale)) scales$y <- yscale
    }
    
    if(!missing(groups)) groups <- eval(substitute(groups),data)
    
    if(!missing(subset)) subset <- eval(substitute(subset),data)
    
    dul <- options('drop.unused.levels')
    options(drop.unused.levels=FALSE)   ## for empty cells
    on.exit(options(dul))                      ## across some panels
    
    do.call("xyplot", c(list(x = formula, data=data, prepanel=prepanel,
                             panel=panel, horizontal = horizontal),
                        if(length(ylab))list(ylab=ylab),
                        if(length(ylim))list(ylim=ylim),
                        if(length(xlab))list(xlab=xlab),
                        if(!missing(groups))list(groups=groups),
                        if(!missing(subset))list(subset=subset),
                        if(length(scales))list(scales=scales),
                        list(...)))
  }
