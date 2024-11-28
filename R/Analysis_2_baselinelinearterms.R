##########################################
### Baseline analysis: effect of T, ah ###
##########################################
Fridaydata <- alldata_weekly_splitJanssen %>%
   filter(dayinweek == "Friday")


regressionobjects_2 <- list()
regressionsummaries_2 <- tibble()

modelid <- 1
for(ah in 1:2) {
  for(temp in 1:2) {
    for(delayedX in 1:2) {
      for(offset in 1:2) {
        vergelijking <- paste("logRt ~ period",
                              c("", "+ offset(logUnprotected + logRvariant)")[offset],
                              c("", "+ ah")[ah],
                              c("", "+ temp")[temp],
                              c("", "", "", "+ ah_1")[ah + delayedX],
                              c("", "", "", "+ temp_1")[temp + delayedX])
        gamresult <- gam(formula = as.formula(vergelijking), data = Fridaydata, weights = weight, method = "REML")
        
        regressionobjects_2 <- c(regressionobjects_2, list(gamresult))
        
        intervals <- simulate_intervals(gamresult)
        
        if(ah == 2) {
          ah_0 <- present_interval(intervals, "ah")
        } else {
          ah_0 <- ""
        }
        if(temp == 2) {
          temp_0 <- present_interval(intervals, "temp")
        } else {
          temp_0 <- ""
        }
        if(ah == 2 && delayedX == 2) {
          ah_1 <- present_interval(intervals, "ah_1")
        } else {
          ah_1 <- ""
        }
        if(temp == 2 && delayedX == 2) {
          temp_1 <- present_interval(intervals, "temp_1")
        } else {
          temp_1 <- ""
        }
        
        npar <- attr(logLik(gamresult),"df")
        nsample <- nobs(gamresult)
        
        regressionsummaries_2 <- bind_rows(
          regressionsummaries_2,
          tibble(modelid = modelid,
                 offset = (offset == 2),
                 Xdelay = (delayedX == 2),
                 ah = (ah == 2),
                 temp = (temp == 2),
                 ah_0 = ah_0,
                 temp_0 = temp_0,
                 ah_1 = ah_1,
                 temp_1 = temp_1,
                 AIC = AIC(gamresult),
                 AICc = AIC(gamresult) + 2 * npar * (npar + 1) / (nsample - npar - 1),
                 BIC = BIC(gamresult),
                 Rsq = R_squared(gamresult))
        )
        modelid <- modelid + 1
      }
    }
  }
}

