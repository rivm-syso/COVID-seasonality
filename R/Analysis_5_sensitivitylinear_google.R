##########################################
### Baseline analysis: effect of T, ah ###
##########################################
Fridaydata <- alldata_weekly_splitJanssen %>%
  filter(dayinweek == "Friday") 


regressionobjects_5 <- list()
regressionsummaries_5 <- tibble()


modelid <- 1
for(climate in c("", "+ temp", "+ temp + ah")) {
  for(retail in 1:2) {
    for(grocery in 1:2) {
      for(parks in 1:2) {
        for(transit in 1:2) {
          for(work in 1:2) {
            for(residential in 1:2) {
              for(delayedX in 1:2) {
                vergelijking <- paste("logRt ~ offset(logUnprotected + logRvariant) + period",
                                      climate,
                                      c("", "+ retail")[retail],
                                      c("", "+ grocery")[grocery],
                                      c("", "+ parks")[parks],
                                      c("", "+ transit")[transit],
                                      c("", "+ work")[work],
                                      c("", "+ residential")[residential],
                                      c("", "", "", "+ retail_1")[retail + delayedX],
                                      c("", "", "", "+ grocery_1")[grocery + delayedX],
                                      c("", "", "", "+ parks_1")[parks + delayedX],
                                      c("", "", "", "+ transit_1")[transit + delayedX],
                                      c("", "", "", "+ work_1")[work + delayedX],
                                      c("", "", "", "+ residential_1")[residential + delayedX])
                gamresult <- gam(formula = as.formula(vergelijking), data = Fridaydata, weights = weight, method = "REML")
                
                regressionobjects_5 <- c(regressionobjects_5, list(gamresult))
                
                npar <- attr(logLik(gamresult),"df")
                nsample <- nobs(gamresult)
                
                regressionsummaries_5 <- bind_rows(
                  regressionsummaries_5,
                  tibble(modelid = modelid,
                         climate = climate, 
                         Xdelay = (delayedX == 2),
                         retail = (retail == 2),
                         grocery = (grocery == 2),
                         parks = (parks == 2),
                         transit = (transit == 2),
                         work = (work == 2),
                         residential = (residential == 2),
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
}

regressionsummaries_5 <- regressionsummaries_5 %>%
  mutate(AkaikeWeight = exp(-AICc/2)/sum(exp(-AICc/2)),
         dAICc = AICc - AICc[1]) 

