# library(splines)
# library(data.table)

mode <- function( x ) {
  ux <- unique( x )
  ux[ which.max( tabulate( match( x, ux ) ) ) ]
}

lrtesticenReg <- function( fit, vars = NULL ) {
  if ( is.null( vars ) )
    vars <- 1:length( attr( terms( fit$formula ), "term.labels" ) )
  res <- do.call( rbind, lapply( vars, function( i ) {
    fit2 <- icenReg::ic_sp( terms( fit$formula )[ -i ], data = fit$getRawData() )
    data.frame( chisq = 2 * ( fit$llk - fit2$llk ), df = length( coef( fit ) ) - length( coef( fit2 ) ),
                p = 1 - pchisq( 2 * ( fit$llk - fit2$llk ), length( coef( fit ) ) - length( coef( fit2 ) ) ) )
  } ) )
  # rownames( res ) <- all.vars( fit$formula[ -1 ] )
  rownames( res ) <- attr( terms( fit$formula )[ vars ], "term.labels" )
  res
}

pred <- function( fit, vars = NULL, conf.int = 0.95, cont.quant.limit = c( 0.05, 0.95 ), cont.len = 1000 ) {
  RawData <- fit$getRawData()
  bs_smps <- if( fit$par=="semi-parametric" ) fit$bsMat else do.call( rbind, fit$mcmcList )[ , -(1:2) ]
  if ( is.null( vars ) ) {
    vars <- all.vars( fit$formula[ -2 ] )
  }
  vardef <- lapply( all.vars( fit$formula[ -2 ] ),
                    function( x ) if( length( unique( x ) )>6 ) median( RawData[[ x ]] ) else mode( RawData[[ x ]] ) )
  names( vardef ) <- all.vars( fit$formula[ -2 ] )
  preds <- lapply( vars, function( var ) {
    predgrid <- vardef
    if( length( unique( RawData[[ var ]] ) )>6 ) {
      predgrid[[ var ]] <- seq( quantile( RawData[[ var ]], cont.quant.limit[ 1 ], na.rm = TRUE ),
                                quantile( RawData[[ var ]], cont.quant.limit[ 2 ], na.rm = TRUE ), length.out = cont.len )
    } else {
      predgrid[[ var ]] <- unique( na.omit( RawData[[ var ]] ) )
    }
    predgrid <- do.call( expand.grid, predgrid )
    predgridcov <- icenReg:::expandX( fit$formula, data = predgrid, fit ) 
    predgridcov <- predgridcov - matrix( rep( fit$covarOffset, nrow( predgridcov ) ), nr = nrow( predgridcov ), byrow = TRUE )
    all_bs_lps <- bs_smps %*% t( predgridcov )
    res <- data.frame( #x = if( is.factor( predgrid[[ var ]] ) ) as.numeric(predgrid[[ var ]]) else predgrid[[ var ]],
                       #xlab = if(is.factor(predgrid[[var]])) as.character(predgrid[[var]]) else rep(NA_character_, length(predgrid[[var]])),
                       x = as.character(predgrid[[ var ]]),
                       var = var,
                       pred = if ( fit$par=="semi-parametric" )
                         matrix( rep( predict( fit, type = "lp", predgrid ), 3 ), nc = 3 ) +
                         sqrt( diag( predgridcov %*% vcov( fit ) %*% t( predgridcov ) ) ) %*% t( c( 0, -1, 1 ) * qnorm( 1-conf.int/2 ) )
                       else cbind( predict( fit, type = "lp", predgrid ), NA, NA ),
                       t( apply( all_bs_lps, 2, function( x ) quantile( x, p = c( 0.5, ( 1-conf.int )/2, 1-( 1-conf.int )/2 ) ) ) ),
                       stringsAsFactors = FALSE )
    names( res )[ 3:8 ] <- c( "pred", "normcilwr", "normciupr", "bspred", "bscilwr", "bsciupr" )
    res
  } )
  preds <- do.call( rbind, preds )
  preds$var <- factor( preds$var, unique( preds$var ) )
  return( preds )
}

plotpred <- function( preds, ci.type = "normal", lrtest = NULL, yfree = FALSE, discrete.hr.plot = TRUE, xlab = "x" ) {
  preds <- if ( ci.type=="normal" ) preds[ , c( "x", "var", "pred", "normcilwr", "normciupr" ) ] else
    preds[ , c( "x", "var", "pred", "bscilwr", "bsciupr" ) ]
  names( preds )[ 4:5 ] <- c( "lwr", "upr" )
  Hmisc::xYplot( Hmisc::Cbind( pred, lwr, upr ) ~ x | var, data = preds, ylab = "log hazard", xlab = xlab,
                 ylim = if( yfree ) lapply( unique( preds$var ), function( v )
                   extendrange( c( preds$lwr[ preds$var==v ], preds$upr[ preds$var==v ] ) ) ) else c( NULL, NULL ),
                 prepanel = function( x, y, ... ) {
                   if( length( unique( x ) )>6 )
                     list( xlim = range( as.numeric(x), na.rm = TRUE ), ylim = extendrange( y ) )
                   else
                     list( xlim = levels( as.factor( x ) ), ylim = extendrange( y ) )
                 },
                 panel = function( x, y, ...) {
                   if( length( unique( x ) )>6 ){
                     Hmisc::panel.xYplot( as.numeric(x), y, type = "l", method = "filled bands", col.fill = gray( seq( .825, .55, length = 5 ) ), ... )
                   }
                   else {
                     Hmisc::panel.Dotplot( as.factor( x ), y, horizontal = FALSE, ... )
                     if( discrete.hr.plot ) lattice::panel.text( as.factor( x )[ -1 ], y[ -1 ], round( ( exp( y-y[ 1 ] ) )[ -1 ], 2 ), pos = 4 )
                   }
                   if( !is.null( lrtest ) ) {
                     pval <- Hmisc::format.pval( lrtest[ lattice::panel.number(), 3 ], digits = 3, eps = 10^(-3))
                     pval <- ifelse( grepl( "<", pval ), paste( "P", pval, sep = "" ), paste( "P==", pval, sep = "" ) )
                     cpl <- lattice::current.panel.limits()
                     lattice::ltext( mean( cpl$xlim ), 0.9*cpl$ylim[2]+0.1*cpl$ylim[1],
                                     parse( text = paste( "chi[",  lrtest[ lattice::panel.number(), 2 ], "]^2 == ",
                                                          round( lrtest[ lattice::panel.number(), 1 ], 1 ), "~~", pval ) ) )
                   }
                 }, scales = list( relation = "free" ) )
}

plotsummary <- function(fit) {
  res <- data.frame(confint(fit), fit$coefficients)
  colnames(res) <- c("lwr", "upr", "HR")
  res <- exp(res)
  vars <- all.vars(fit$formula)[-(1:2)]
  res$varlabs <- rep(vars, sapply(vars, function(x) if(is.factor(fit$.dataEnv$data[[x]])) length(levels(fit$.dataEnv$data[[x]]))-1 else 1))
  res$reflabs <- as.vector(unlist(sapply(vars, function(x) if(is.factor(fit$.dataEnv$data[[x]]))
    paste0(levels(fit$.dataEnv$data[[x]])[-1], ":", levels(fit$.dataEnv$data[[x]])[1]) else "+1")))
  res$lab <- paste0(res$varlabs, " - ", res$reflabs)
  res$lab <- factor(res$lab, levels = res$lab)
  ggplot2::ggplot(res, ggplot2::aes(x = HR, xmin = lwr, xmax = upr, y = lab)) + ggplot2::geom_point() +
    ggplot2::geom_errorbar(width = 0.5) + ggplot2::labs(y = "") + ggplot2::scale_x_log10() + ggplot2::coord_cartesian(xlim = c(0.1, 500)) +
    ggplot2::geom_vline(xintercept = 1, color = "red") + ggplot2::scale_y_discrete(limits = rev(levels(res$lab)))
}