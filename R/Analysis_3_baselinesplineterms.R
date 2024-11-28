##########################################
### Baseline analysis: effect of T, ah ###
##########################################
Fridaydata <- alldata_weekly_splitJanssen %>%
  filter(dayinweek == "Friday")  


regressionobjects_3 <- list()
regressionsummaries_3 <- tibble()

modelid <- 1
for(ah in 1:2) {
  for(temp in 1:2) {
    for(delayedX in 1:2) {
      for(offset in 1:2) {
        vergelijking <- paste("logRt ~ period",
                              c("", "+ offset(logUnprotected + logRvariant)")[offset],
                              c("", "+ s(ah)")[ah],
                              c("", "+ s(temp)")[temp],
                              c("", "", "", "+ s(ah_1)")[ah + delayedX],
                              c("", "", "", "+ s(temp_1)")[temp + delayedX])
        gamresult <- gam(formula = as.formula(vergelijking), data = Fridaydata, weights = weight, method = "REML")
        
        regressionobjects_3 <- c(regressionobjects_3, list(gamresult))
        
        intervals <- simulate_intervals(gamresult)
        
        if(ah == 2) {
          ah_0 <- pen.edf(gamresult)["s(ah)"] 
        } else {
          ah_0 <- 0
        }
        if(temp == 2) {
          temp_0 <- pen.edf(gamresult)["s(temp)"] 
        } else {
          temp_0 <- 0
        }
        if(ah == 2 && delayedX == 2) {
          ah_1 <- pen.edf(gamresult)["s(ah_1)"] 
        } else {
          ah_1 <- 0
        }
        if(temp == 2 && delayedX == 2) {
          temp_1 <- pen.edf(gamresult)["s(temp_1)"]
        } else {
          temp_1 <- 0
        }
        
        npar <- attr(logLik(gamresult),"df")
        nsample <- nobs(gamresult)
        
        regressionsummaries_3 <- bind_rows(
          regressionsummaries_3,
          tibble(modelid = modelid,
                 offset = (offset == 2),
                 Xdelay = (delayedX == 2),
                 ah = (ah == 2),
                 temp = (temp == 2),
                 ah_0 = round(ah_0, 2),
                 temp_0 = round(temp_0, 2),
                 ah_1 = round(ah_1, 2),
                 temp_1 = round(temp_1, 2),
                 npar = npar,
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


