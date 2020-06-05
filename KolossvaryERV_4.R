library(icenReg)
library(doParallel)
library(splines)
library(data.table)

source("plotpredicticenRegRoutines.R")

RawData <- foreign::read.dta("ERVProspektiv20180121corrected.dta")
cbind( RawData$Betegaz[is.na(RawData$meghalt)])
RawData <- RawData[!is.na(RawData$meghalt),]
RawData$Betegaz[is.na( RawData$DATEV1)]
Reduce(union, list(RawData$Betegaz[!is.na(RawData$DATEV1) & !is.na(RawData$DATEV2) & (RawData$DATEV1 > RawData$DATEV2)],
                   RawData$Betegaz[!is.na(RawData$DATEV2) & !is.na(RawData$DATEV3) & (RawData$DATEV2 > RawData$DATEV3)],
                   RawData$Betegaz[!is.na(RawData$DATEV3) & !is.na(RawData$DATEV4) & (RawData$DATEV3 > RawData$DATEV4)],
                   RawData$Betegaz[!is.na(RawData$DATEV4) & !is.na(RawData$DATEV5) & (RawData$DATEV4 > RawData$DATEV5)],
                   RawData$Betegaz[!is.na(RawData$DATEV5) & !is.na(RawData$DATEV6) & (RawData$DATEV5 > RawData$DATEV6)]))

prop.table(table(RawData$ScoreV1, RawData$PAD09V1), 1)[, 2]*100
prop.table(table(RawData$ScoreV1==4, RawData$PAD09V1), 1)[, 2]*100

RawData$LastVis <- as.numeric(apply(RawData[, c("DATEV1", "DATEV2", "DATEV3", "DATEV4", "DATEV5", "DATEV6")], 1,
                                    function(x) tail(which(!is.na(x)), 1)))
RawData$LastVisDate <- as.Date(as.character(apply(RawData[, c("DATEV1", "DATEV2", "DATEV3", "DATEV4", "DATEV5", "DATEV6")], 1,
                                                  function(x) tail(x[!is.na(x)], 1)), origin = "1970-01-01"))
RawData$BeforeLastVisDate <- as.Date(as.character(apply(RawData[, c("DATEV1", "DATEV2", "DATEV3", "DATEV4", "DATEV5", "DATEV6")],
                                                        1, function(x) tail(x[!is.na(x)], 2)[1]), origin = "1970-01-01"))

RawData$Right <- as.Date(NA)
RawData$Left <- RawData$LastVisDate
RawData$Right[RawData$meghalt==1] <- RawData$LastVisDate[RawData$meghalt==1]
RawData$Left[RawData$meghalt==1] <- RawData$BeforeLastVisDate[ RawData$meghalt==1]

RawData$Right <- as.numeric(RawData$Right-RawData$DATEV1)
RawData$Left <- as.numeric(RawData$Left-RawData$DATEV1)

RawData$Right2 <- as.Date(NA)
RawData$Right2[RawData$meghalt==1] <- RawData$LastVisDate[RawData$meghalt==1]
RawData$Left2 <- RawData$LastVisDate

RawData$Right2 <- as.numeric(RawData$Right2-RawData$DATEV1)
RawData$Left2 <- as.numeric(RawData$Left2-RawData$DATEV1)

for(var in c("Neme", "PADgroupV1", "PAD09V1", "ScoreV1", "SszelutV1", "SerszV1", "SinfarktV1", "DIABCOMPLV1", "GFRV1",
             "GFR60V1")) {
  RawData[[var]] <- as.factor(RawData[[var]])
}

RawData <- RawData[, c("Left", "Right", "Neme", "AGE", "DIABCOMPLV1", "SinfarktV1", "SszelutV1",
                         "BMIV1", "GFR60V1", "PAD09V1", "ScoreV1", "RRsystV1", "RRdiastV1")]
RawData <- as.data.table(RawData)
RawData <- na.omit(RawData, cols = c(1, 3:ncol(RawData)))

myClust <- makeCluster(detectCores()-1)
registerDoParallel(myClust)
fit1 <- ic_sp(Surv(Left, Right, type = "interval2") ~ Neme + ns(AGE, df = 3) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                 ns(BMIV1, df = 3) + GFR60V1 + PAD09V1 + ns(RRsystV1, df = 3), data = RawData,
               model = "ph", bs_samples = 1000, useMCores = TRUE)
fit2A <- ic_sp( Surv(Left, Right, type = "interval2") ~ Neme + ns(AGE, df = 3) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                  ns(BMIV1, df = 3) + GFR60V1 + PAD09V1 + ns(RRsystV1, df = 3) + ns( RRdiastV1, df = 3), data = RawData,
                model = "ph", bs_samples = 1000, useMCores = TRUE)
fit2B <- ic_sp( Surv(Left, Right, type = "interval2") ~ Neme + ns(AGE, df = 3) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                  ns(BMIV1, df = 3) + GFR60V1 + PAD09V1 * ns(RRsystV1, df = 3), data = RawData,
                model = "ph", bs_samples = 1000, useMCores = TRUE)
fit2C <- ic_sp( Surv(Left, Right, type = "interval2") ~ Neme + ns(AGE, df = 3) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                  ns(BMIV1, df = 3) + GFR60V1 + PAD09V1 + ns(RRsystV1, df = 3) * ns(RRdiastV1, df = 3), data = RawData,
                model = "ph", bs_samples = 1000, useMCores = TRUE)
fit3 <- ic_sp( Surv(Left, Right, type = "interval2") ~ ScoreV1 + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                 ns(BMIV1, df = 3) + GFR60V1 + PAD09V1, data = RawData,
               model = "ph", bs_samples = 1000, useMCores = TRUE)
fit4 <- ic_sp( Surv(Left, Right, type = "interval2") ~ ns(AGE, df = 3) + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                 ns(BMIV1, df = 3) + GFR60V1 + Neme * PAD09V1 + ns(RRsystV1, df = 3), data = RawData,
               model = "ph", bs_samples = 1000, useMCores = TRUE)
fit5 <- ic_sp( Surv(Left, Right, type = "interval2") ~ Neme + DIABCOMPLV1 + SinfarktV1 + SszelutV1 +
                 ns(BMIV1, df = 3) + GFR60V1 + PAD09V1 * ns(AGE, df = 3) + ns(RRsystV1, df = 3), data = RawData,
               model = "ph", bs_samples = 1000, useMCores = TRUE)
fit6 <- ic_sp( Surv(Left, Right, type = "interval2") ~ ScoreV1, data = RawData,
               model = "ph", bs_samples = 1000, useMCores = TRUE)
fit6B <- ic_sp( Surv(Left, Right, type = "interval2") ~ ScoreV1 * PAD09V1, data = RawData,
                model = "ph", bs_samples = 1000, useMCores = TRUE)
fit6C <- ic_sp( Surv(Left, Right, type = "interval2") ~ ScoreV1 + PAD09V1, data = RawData,
                model = "ph", bs_samples = 1000, useMCores = TRUE)
stopCluster(myClust)

#models <- c( fit1 = fit1, fit2A = fit2A, fit2B = fit2B, fit2C = fit2C, fit3 = fit3 )
models <- c( fit4 = fit4, fit5 = fit5, fit6 = fit6, fit6B = fit6B, fit6C = fit6C )

for(i in 1:length(models)) {
  temp <- data.frame(HR = exp(coef(models[[i]])), exp(confint(models[[i]])),
                      p = Hmisc::format.pval(summary(models[[i]])$summaryParameters[, 5], eps = 0.001))
  write.csv2(temp, paste0(names(models)[i], ".csv"))
}

predfit <- pred(fit1)
lrtestfit <- lrtesticenReg(fit1)
plotpred(predfit, ci.type = "normal",lrtest = lrtestfit)
plotpred(predfit, ci.type = "bs",lrtest = lrtestfit)
