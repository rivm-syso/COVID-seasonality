###############################################
### pre- Baseline analysis: autocorrelation ###
###############################################
Fridaydata_mergedJanssen <- alldata_weekly %>% 
  filter(dayinweek == "Friday")
Fridaydata <- alldata_weekly_splitJanssen %>% 
  filter(dayinweek == "Friday")

regressionobjects_1 <- list()
regressionsummaries_1 <- tibble()

modelid <- 1
for(daily in 1:2) {
  for(delay in 1:2) {
    for(season in 1:2) {
      vergelijking <- paste("logRt ~ offset(logUnprotected) + offset(logRvariant) + period",
                            c("", "+ dayinweek")[daily],
                            c("", "", " + logRt_1")[delay + 2 - daily],
                            c("", "", "", paste0(" + logRt_", 1:7, collapse = ""))[delay + daily],
                            c("", "+ ah + temp")[season])
      if(daily == 2) {
        lmresult <- lm(formula = vergelijking, data = alldata_daily_splitJanssen, weights = weight)
      } else {
        lmresult <- lm(formula = vergelijking, data = Fridaydata, weights = weight)
        
      }
      regressionobjects_1 <- c(regressionobjects_1, list(lmresult))
      correlationtest <- cor.test(head(lmresult$residuals, -1), tail(lmresult$residuals, -1))
      regressionsummaries_1 <- bind_rows(
        regressionsummaries_1,
        tibble(modelid = modelid,
               DelayTerms = (delay == 2),
               DailyData = (daily == 2),
               Seasonality = (season == 2),
               AIC = AIC(lmresult),
               BreuschGodfrey_p = round(lmtest::bgtest(lmresult)$p.value, 3))
      )
      modelid <- modelid + 1
    }
  }
}



# Conclusion: weekly data, no autocorrelation term

#####################################
### pre-Baseline: detect outliers ###
#####################################

include4cook <- Fridaydata_mergedJanssen$period %in% 
  unique(Fridaydata_mergedJanssen$period[duplicated(Fridaydata_mergedJanssen$period)])

cooksplotdata <- tibble()

for(modelid in 1:2) {
  cooksdistances <- cooks.distance(regressionobjects_1[[modelid]])
  cooksdistances[!include4cook] <- NA
  dataset <- tibble(
    modelid = modelid,
    date = Fridaydata_mergedJanssen$date,
    cooksdist = cooksdistances,
    logRt = Fridaydata_mergedJanssen$logRt,
    temp = Fridaydata_mergedJanssen$temp,
    ah = Fridaydata_mergedJanssen$ah
  )
  cooksplotdata <- bind_rows(cooksplotdata, dataset)
}

include4cook <- Fridaydata_mergedJanssen$period %in% 
  unique(Fridaydata_mergedJanssen$period[duplicated(Fridaydata_mergedJanssen$period)])

for(modelid in 3:4) {
  vergelijking <- paste("logRt ~ offset(logUnprotected) + offset(logRvariant) + period",
                        c("", "+ ah + temp")[modelid - 2])
  
  lmresult <- lm(formula = vergelijking, data = Fridaydata_mergedJanssen, weights = weight)

  cooksdistances <- cooks.distance(lmresult)
  cooksdistances[!include4cook] <- NA
  dataset <- tibble(
    modelid = modelid,
    date = Fridaydata_mergedJanssen$date,
    cooksdist = cooksdistances,
    logRt = Fridaydata_mergedJanssen$logRt,
    temp = Fridaydata_mergedJanssen$temp,
    ah = Fridaydata_mergedJanssen$ah
  )
  cooksplotdata <- bind_rows(cooksplotdata, dataset)
}

# Conclusion: split period 26 into two separate weeks

