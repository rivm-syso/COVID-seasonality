##########################################
### Baseline analysis: effect of T, ah ###
##########################################
Fridaydata <- alldata_weekly_splitJanssen %>%
   filter(dayinweek == "Friday")


regressionobjects_10 <- list()
regressionsummaries_10 <- tibble()

modelid <- 1
for(ah in 1:2) {
  for(temp in 1:2) {
    for(delayedX in 1:2) {
      for(dayofweek in c("Friday", "Saturday", "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday")) {
        vergelijking <- paste("logRt ~ period + offset(logUnprotected + logRvariant)",
                              c("", "+ s(ah)")[ah],
                              c("", "+ s(temp)")[temp],
                              c("", "", "", "+ s(ah_1)")[ah + delayedX],
                              c("", "", "", "+ s(temp_1)")[temp + delayedX])
        gamresult <- gam(formula = as.formula(vergelijking), data = alldata_weekly_splitJanssen %>%
                           filter(dayinweek == dayofweek), weights = weight, method = "REML")
        
        regressionobjects_10 <- c(regressionobjects_10, list(gamresult))
        
        relRt <- summarise_relRt_season(gamresult)
        relRt <- present_interval(t(relRt), "amplitude", precision = 1)

        npar <- attr(logLik(gamresult),"df")
        nsample <- nobs(gamresult)
        
        regressionsummaries_10 <- bind_rows(
          regressionsummaries_10,
          tibble(modelid = modelid,
                 day = dayofweek,
                 Xdelay = (delayedX == 2),
                 ah = (ah == 2),
                 temp = (temp == 2),
                 AIC = AIC(gamresult),
                 AICc = AIC(gamresult) + 2 * npar * (npar + 1) / (nsample - npar - 1),
                 BIC = BIC(gamresult),
                 Rsq = R_squared(gamresult),
                 relRt = relRt)
        )
        modelid <- modelid + 1
      }
    }
  }
}

regressionsummaries_10 <- regressionsummaries_10 %>%
  mutate(AkaikeWeight = exp(-AICc/2)/sum(exp(-AICc/2)),
         dAICc = AICc - AICc[1], .by = day) 

regressionsummaries_10 %>%
  filter(!Xdelay | ah | temp) %>%
  dplyr::select(temp, ah, Xdelay, day, dAICc) %>%
  arrange(Xdelay, ah, temp) %>%
  pivot_wider(names_from = day, values_from = dAICc)

regressionsummaries_10 %>%
  filter(!Xdelay | ah | temp) %>%
  dplyr::select(temp, ah, Xdelay, day, relRt) %>%
  arrange(Xdelay, ah, temp) %>%
  pivot_wider(names_from = day, values_from = relRt)
