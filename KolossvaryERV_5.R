library(data.table)
library(splines)

RawData <- haven::read_dta("ERVProspektiv20180121corrected.dta", .name_repair = "universal")
cbind( RawData$Betegaz[is.na(RawData$meghalt)])
RawData <- RawData[!is.na(RawData$meghalt), ]
RawData$Betegaz[is.na( RawData$DATEV1)]
Reduce(union, list(RawData$Betegaz[!is.na(RawData$DATEV1) & !is.na(RawData$DATEV2) & (RawData$DATEV1 > RawData$DATEV2)],
                   RawData$Betegaz[!is.na(RawData$DATEV2) & !is.na(RawData$DATEV3) & (RawData$DATEV2 > RawData$DATEV3)],
                   RawData$Betegaz[!is.na(RawData$DATEV3) & !is.na(RawData$DATEV4) & (RawData$DATEV3 > RawData$DATEV4)],
                   RawData$Betegaz[!is.na(RawData$DATEV4) & !is.na(RawData$DATEV5) & (RawData$DATEV4 > RawData$DATEV5)],
                   RawData$Betegaz[!is.na(RawData$DATEV5) & !is.na(RawData$DATEV6) & (RawData$DATEV5 > RawData$DATEV6)]))

prop.table( table( RawData$ScoreV1, RawData$PAD09V1 ), 1 )[ , 2 ]*100
prop.table( table( RawData$ScoreV1==4, RawData$PAD09V1 ), 1 )[ , 2 ]*100

RawData$LastVis <- as.numeric( apply( RawData[ , c( "DATEV1", "DATEV2", "DATEV3", "DATEV4", "DATEV5", "DATEV6" ) ], 1,
                                      function( x ) tail( which( !is.na( x ) ), 1 ) ) )
RawData$LastVisDate <- as.Date( as.character( apply( RawData[ , c( "DATEV1", "DATEV2", "DATEV3", "DATEV4", "DATEV5", "DATEV6" ) ], 1,
                                                     function( x ) tail( x[ !is.na( x ) ], 1 ) ), origin = "1970-01-01" ) )
RawData$BeforeLastVisDate <- as.Date( as.character( apply( RawData[ , c( "DATEV1", "DATEV2", "DATEV3", "DATEV4", "DATEV5", "DATEV6" ) ], 1,
                                                           function( x ) tail( x[ !is.na( x ) ], 2 )[1] ), origin = "1970-01-01" ) )

RawData$Right <- as.Date( NA )
RawData$Left <- RawData$LastVisDate
RawData$Right[ RawData$meghalt==1 ] <- RawData$LastVisDate[ RawData$meghalt==1 ]
RawData$Left[ RawData$meghalt==1 ] <- RawData$BeforeLastVisDate[ RawData$meghalt==1 ]

RawData$Right <- as.numeric( RawData$Right-RawData$DATEV1 )
RawData$Left <- as.numeric( RawData$Left-RawData$DATEV1 )

RawData$Right2 <- as.Date( NA )
RawData$Right2[ RawData$meghalt==1 ] <- RawData$LastVisDate[ RawData$meghalt==1 ]
RawData$Left2 <- RawData$LastVisDate

RawData$Right2 <- as.numeric( RawData$Right2-RawData$DATEV1 )
RawData$Left2 <- as.numeric( RawData$Left2-RawData$DATEV1 )

for(var in c("Neme", "PADgroupV1", "PAD09V1", "ScoreV1", "SszelutV1", "SerszV1", "SinfarktV1", "DIABCOMPLV1", "GFRV1",
             "GFR60V1")) RawData[[ var ]] <- as.factor( RawData[[ var ]] )

RawData$ABImaxV1 <- (RawData$BABImaxV1 + RawData$JABImaxV1)/2
RawData <- data.table( RawData )

ggplot2::ggplot(melt(RawData[ABImaxV1<2.5, .(Betegaz, `Left side` = BABImaxV1, `Right side` = JABImaxV1)], id.vars = "Betegaz"),
                ggplot2::aes(x = value, color = variable, group = variable)) + ggplot2::geom_density(adjust = 2) +
  ggplot2::labs(x = "Ankle brachial index", y = "Density", color = "")
ggplot2::ggsave("ABIEloszlasKDE.pdf", width = 16, height = 9)

RawData <- RawData[ , c( "Left", "Right", "Neme", "AGE", "DIABCOMPLV1", "SinfarktV1", "SszelutV1",
                         "BMIV1", "GFR60V1", "PAD09V1", "ScoreV1", "RRsystV1", "RRdiastV1", "ABImaxV1") ]

RawData <- na.omit(RawData, cols = c(1, 3:ncol(RawData)))

myClust <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(myClust)

fit1 <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ Neme + ns( AGE, df = 3 ) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                         #ns( BMIV1, df = 3 ) + GFR60V1 + PAD09V1 + ns( RRsystV1, df = 3 ), data = RawData,
                         ns( BMIV1, df = 3 ) + GFR60V1 + ns(ABImaxV1, df = 3) + ns( RRsystV1, df = 3 ), data = RawData,
                       model = "ph", bs_samples = 1000, useMCores = TRUE)
fit2A <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ Neme + ns( AGE, df = 3 ) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                          ns( BMIV1, df = 3 ) + GFR60V1 + ns(ABImaxV1, df = 3) + ns( RRsystV1, df = 3 ) + ns( RRdiastV1, df = 3 ), data = RawData,
                        model = "ph", bs_samples = 1000, useMCores = TRUE)
fit2AInt <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ ns( AGE, df = 3 ) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                             ns( BMIV1, df = 3 ) + GFR60V1 + Neme * ns(ABImaxV1, df = 3) + ns( RRsystV1, df = 3 ) + ns( RRdiastV1, df = 3 ), data = RawData,
                           model = "ph", bs_samples = 1000, useMCores = TRUE)
fit2B <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ Neme + ns( AGE, df = 3 ) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                          ns( BMIV1, df = 3 ) + GFR60V1 + ns(ABImaxV1, df = 3) * ns( RRsystV1, df = 3 ), data = RawData,
                        model = "ph", bs_samples = 1000, useMCores = TRUE)
fit2C <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ Neme + ns( AGE, df = 3 ) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                          ns( BMIV1, df = 3 ) + GFR60V1 + ns(ABImaxV1, df = 3) + ns( RRsystV1, df = 3 ) * ns( RRdiastV1, df = 3 ), data = RawData,
                        model = "ph", bs_samples = 1000, useMCores = TRUE)
fit3 <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ ScoreV1 + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                         ns( BMIV1, df = 3 ) + GFR60V1 + ns(ABImaxV1, df = 3), data = RawData,
                       model = "ph", bs_samples = 1000, useMCores = TRUE)
fit4 <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ ns( AGE, df = 3 ) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                         ns( BMIV1, df = 3 ) + GFR60V1 + Neme * ns(ABImaxV1, df = 3) + ns( RRsystV1, df = 3 ), data = RawData,
                       model = "ph", bs_samples = 1000, useMCores = TRUE)
fit5 <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ Neme + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                         ns( BMIV1, df = 3 ) + GFR60V1 + ns(ABImaxV1, df = 3) * ns( AGE, df = 3 ) + ns( RRsystV1, df = 3 ), data = RawData,
                       model = "ph", bs_samples = 1000, useMCores = TRUE)
fit6 <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ ScoreV1, data = RawData,
                       model = "ph", bs_samples = 1000, useMCores = TRUE)
fit6B <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ ScoreV1 * ns(ABImaxV1, df = 3), data = RawData,
                        model = "ph", bs_samples = 1000, useMCores = TRUE)
fit6C <- icenReg::ic_sp(survival::Surv(Left, Right, type = "interval2") ~ ScoreV1 + ns(ABImaxV1, df = 3), data = RawData,
                        model = "ph", bs_samples = 1000, useMCores = TRUE)

parallel::stopCluster(myClust)

fits <- list(fit1, fit2A, fit2B, fit2C, fit3, fit4, fit5, fit6, fit6B, fit6C)

saveRDS(fits, "ERVfits.rds")

# fits <- readRDS("ERVfits.rds")

for( i in 1:length( fits )  ) {
  temp <- data.frame( HR = exp( coef( fits[[ i ]] ) ), exp( confint( fits[[ i ]] ) ),
                      p = preport(  summary( fits[[ i ]] )$summaryParameters[ , 5 ] ) )
  temp <- temp[ !grepl( "ns(", row.names( temp ), fixed = TRUE ), ]
  temp$Variable <- row.names(temp)
  fwrite(temp, paste0("Modell_20220827_", i, ".csv"))
}

source("plotpredicticenReg.R")

cl <- parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterExport(cl, c("pred", "plotsummary", "plotpred", "lrtesticenReg", "mode", "fits"))
parallel::clusterEvalQ(cl, library(splines))
parallel::clusterEvalQ(cl, library(data.table))

parallel::parLapply(cl, 1:length(fits), function(i) {
  predfit <- pred(fits[[i]], cont.quant.limit = c(0.01, 0.99))
  lrtestfit <- lrtesticenReg(fits[[i]])
  cairo_pdf(paste0("Modell_20220827_", i, ".pdf"), onefile = TRUE, width = 16, height = 9)
  print(plotpred(predfit, ci.type = "bs", lrtest = lrtestfit))
  print(plotpred(predfit, ci.type = "bs", lrtest = lrtestfit, yfree = TRUE))
  dev.off()
})

parallel::stopCluster(cl)