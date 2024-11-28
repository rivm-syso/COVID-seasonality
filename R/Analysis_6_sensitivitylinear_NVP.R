##########################################
### Baseline analysis: effect of T, ah ###
##########################################
Fridaydata <- alldata_weekly_splitJanssen %>%
  filter(dayinweek == "Friday") 


regressionobjects_6 <- list()
regressionsummaries_6 <- tibble()

modelid <- 1
for(climate in c("", "+ temp", "+ temp + ah")) {
  for(thuis in 1:2) {
    for(rondje in 1:2) {
      for(supermarkt in 1:2) {
        for(kantoor in 1:2) {
          for(zorg in 1:2) {
            for(delayedX in 1:2) {
              vergelijking <- paste("logRt ~ offset(logUnprotected + logRvariant) + period",
                                    climate,
                                    c("", "+ index_thuisblijvers")[thuis],
                                    c("", "+ index_rondje")[rondje],
                                    c("", "+ index_supermarkt")[supermarkt],
                                    c("", "+ index_kantoor")[kantoor],
                                    c("", "+ index_zorg")[zorg],
                                    c("", "", "", "+ index_thuisblijvers_1")[thuis + delayedX],
                                    c("", "", "", "+ index_rondje_1")[rondje + delayedX],
                                    c("", "", "", "+ index_supermarkt_1")[supermarkt + delayedX],
                                    c("", "", "", "+ index_kantoor_1")[kantoor + delayedX],
                                    c("", "", "", "+ index_zorg_1")[zorg + delayedX])
              gamresult <- gam(formula = as.formula(vergelijking), data = Fridaydata, weights = weight, method = "REML")
              
              regressionobjects_6 <- c(regressionobjects_6, list(gamresult))
              
              npar <- attr(logLik(gamresult),"df")
              nsample <- nobs(gamresult)
              
              regressionsummaries_6 <- bind_rows(
                regressionsummaries_6,
                tibble(modelid = modelid,
                       climate = climate, 
                       Xdelay = (delayedX == 2),
                       thuis = (thuis == 2),
                       rondje = (rondje == 2),
                       supermarkt = (supermarkt == 2),
                       kantoor = (kantoor == 2),
                       zorg = (zorg == 2),
                       # ah_0 = ah_0,
                       # temp_0 = temp_0,
                       AIC = AIC(gamresult),
                       AICc = AIC(gamresult) + 2 * npar * (npar + 1) / (nsample - npar - 1),
                       BIC = BIC(gamresult),
                       npar = npar,
                       nsample = nsample,
                       Rsq = R_squared(gamresult))
              )
              modelid <- modelid + 1
            }
          }
        }
      }
    }
  }
}



regressionsummaries_6 <- regressionsummaries_6 %>%
  mutate(AkaikeWeight = exp(-AICc/2)/sum(exp(-AICc/2)),
         dAICc = AICc - AICc[1])  

