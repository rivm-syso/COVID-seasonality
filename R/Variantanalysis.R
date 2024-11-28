###################################
### read and clean variant data ###
###################################
Variantdata_raw <- read_csv2("data/COVID-19_varianten_20240814.csv")

Variantdata_clean <- Variantdata_raw %>%
  ### up to Delta
  filter(Date_of_statistics_week_start <= ymd("2021-10-05")) %>% 

  ### remove records of variants arising later
  group_by(Variant_name) %>%
  filter(sum(Variant_cases) > 0) %>%
  ungroup() %>%
  
  ### count number of samples without variant name (most are wild type)
  ### and remove spurious wild type samples after dominance of Alpha
  dplyr::select(Date_of_statistics_week_start, Variant_name, Variant_cases, Sample_size) %>%
  group_by(Date_of_statistics_week_start) %>%
  mutate(WildType = Sample_size - sum(Variant_cases)) %>%
  ungroup() %>%
  mutate(WildType = if_else(Date_of_statistics_week_start > ymd("2021-05-01"), 0, WildType)) %>%
  
  ### rearrange data
  pivot_wider(names_from = Variant_name, values_from = Variant_cases) %>%
  mutate(Days_since_start = as.numeric(Date_of_statistics_week_start - ymd("2020-02-12"))) %>%
  dplyr::select(Days_since_start, WildType:Alpha) %>%
  pivot_longer(WildType:Alpha, names_to = "Variant", values_to = "Count") %>%
  mutate(
    Variant = factor(Variant, levels = c("WildType", "Alpha", "Beta", "Gamma", "Delta"))
    ) 
  
#############################
### fit multinomial model ###
#############################
library(nnet)
variantmodelfit <- multinom(Variant ~ Days_since_start, Variantdata_clean, weights = Count)
variantcoefficients <- rbind(c(0, 0), coef(variantmodelfit))
rownames(variantcoefficients)[1] <- "WildType"

########################################################
### build relative infectivity curve due to variants ###
########################################################
relativeinfectivities_variant <-
  ((variantcoefficients[, 2] + 0.875)^4 / (2 * 0.875^3 + variantcoefficients[, 2] * 0.875^2)) /
  ((variantcoefficients[1, 2] + 0.875)^4 / (2 * 0.875^3 + variantcoefficients[1, 2] * 0.875^2))
  
variantproportions <-
  tibble(Days_since_start = rep(1:1262, 5),
         Variant = rep(c("WildType", "Alpha", "Beta", "Gamma", "Delta"), each = 1262)) %>%

  ## reported day of incidence is 3.5 days too early (date is first day of week), 
  ## and 7 days too late (5 days incubation + 2 days to sample), 
  ## therefore curve is shifted 3.5 days to the past
  mutate(FittedRelativeLogIncidence = variantcoefficients[Variant, "(Intercept)"] + 
           (Days_since_start + 3.5) * variantcoefficients[Variant, "Days_since_start"]) %>%
  group_by(Days_since_start) %>%
  
  ## prevent Infinity when taking exponential
  mutate(FittedRelativeLogIncidence = FittedRelativeLogIncidence - max(FittedRelativeLogIncidence),
         
         ## incidence per variant relative to other variants
         FittedRelativeIncidence = exp(FittedRelativeLogIncidence),
         
         ## prevalence per variant relative to other variants (scaled because of infectivity)
         FittedRelativePrevalence = FittedRelativeIncidence / relativeinfectivities_variant[Variant],
         
         ## proportional prevalence per day
         FittedPrevalenceProportion = FittedRelativePrevalence / sum(FittedRelativePrevalence)) %>%
  ungroup() %>%
  dplyr::select(Days_since_start, Variant, FittedPrevalenceProportion)


variants <- variantproportions %>%
  mutate(infectivityweight = relativeinfectivities_variant[Variant]) %>%
  group_by(Days_since_start) %>%
  
  ## sum infectivities of variants, weighed by their prevalences
  summarise(variantcurve = weighted.mean(infectivityweight, FittedPrevalenceProportion)) %>%
  mutate(date = Days_since_start + ymd("2020-02-12")) %>%
  ungroup() %>%
  dplyr::select(date, variantcurve)
  
  
  
