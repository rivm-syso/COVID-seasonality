##########################################
### Baseline analysis: effect of T, ah ###
##########################################
Fridaydata <- alldata_weekly_splitJanssen %>%
  filter(dayinweek == "Friday")  


regressionobjects_4 <- list()
regressionsummaries_4 <- tibble()

modelid <- 1
for(ahtemp in c("_avg","","_04","_08","_12","_16","_20")) {
  Fitdata <- Fridaydata %>%
    mutate(ah = !!as.name(paste0("ah", ahtemp)),
           temp = !!as.name(paste0("temp", ahtemp)))
  for(ah in 1:2) {
    for(temp in 1:2) {
      for(spline in 0:0) {
        vergelijking <- paste("logRt ~ offset(logUnprotected + logRvariant) + period",
                              c("", "+ ah", "", "+ s(ah)")[ah + 2 * spline],
                              c("", "+ temp", "", "+ s(temp)")[temp + 2 * spline])
        gamresult <- gam(formula = as.formula(vergelijking), data = Fitdata, weights = weight, method = "REML")
        
        regressionobjects_4 <- c(regressionobjects_4, list(gamresult))
        
        intervals <- simulate_intervals(gamresult)
        
        if(ah == 2) {
          if(spline == 0) {
          ah_0 <- 1
          } else {
          ah_0 <- pen.edf(gamresult)["s(ah)"]
        }
        } else {
          ah_0 <- 0
        }
        if(temp == 2) { 
          if(spline == 0) {
          temp_0 <- 1
          } else {
            temp_0 <- pen.edf(gamresult)["s(temp)"]
          }
          } else {
          temp_0 <- 0
        }
        
        npar <- attr(logLik(gamresult),"df")
        nsample <- nobs(gamresult)
        
        regressionsummaries_4 <- bind_rows(
          regressionsummaries_4,
          tibble(modelid = modelid,
                 ahtemp = ahtemp,
                 ah = (ah == 2),
                 temp = (temp == 2),
                 ah_0 = ah_0,
                 temp_0 = temp_0,
                 AIC = AIC(gamresult),
                 AICc = AIC(gamresult) + 2 * npar * (npar + 1) / (nsample - npar - 1),
                 BIC = BIC(gamresult))
        )
        modelid <- modelid + 1
      }
    }
  }
}

