#############################
### Read and prepare datasets one by one
#############################
popsize <- 17282160 # population 2019 


#######
#### Rt
#######
# import epicurve reports and epicurve_ZKH of October 30, with instantaneous R and case R (with logvar) 
epicurve_reports <- readRDS(file = paste0("data/epicurve_reports_2022-04-04.rds"))
epicurve_ZKH <- readRDS(file = paste0("data/epicurve_", as.Date("2020-10-30"), ".rds"))

# transition from ZKH to reports on June 12 and add day of week
epicurve <- bind_rows("ziekenhuis" = epicurve_ZKH %>% filter(dates <= as.Date("2020-06-12")),
                      "meldingen" = epicurve_reports %>% rename(dates = date) %>% filter(dates > as.Date("2020-06-12")),
                      .id = "source") %>% 
  mutate(dayinweek = as.factor(weekdays(dates))) %>%
  dplyr::select(dates, instR, instRlogvar, dayinweek) %>%
  rename(date = dates)
rm(epicurve_reports, epicurve_ZKH)

########################
### Cumulative incidence
########################
# susceptibles from NICE hospital data and Pienter serology
susceptibles <- readRDS(paste0("data/prevalence_NICEZKH_2021-07-03.rds")) %>% 
  replace_na(replace = list(incidence_mean = 0)) %>% 
  mutate(susceptibles = popsize - cumsum(incidence_mean)) %>% 
  dplyr::select(date, susceptibles)


############
### variants
############
source("R/Variantanalysis.R")

###############
### vaccination
###############
## VEdelta = 80%-71%: https://www.sciencedirect.com/science/article/pii/S1201971223005349
vaccdata <- read_csv2("data/COVID-19_vaccinatiegraad_per_gemeente_per_week_leeftijd_20230710.csv")

vaccdata <- vaccdata %>%
  filter(Region_level == "Veiligheidsregio") %>% 
  mutate(covprimcomplete = case_when(Coverage_primary_completed == "<=5" ~ "0",
                                     Coverage_primary_completed == ">=95" ~ "100",
                                     Coverage_primary_completed == "9999" ~ NA_character_,
                                     TRUE ~ Coverage_primary_completed),
         covprimcomplete = as.numeric(covprimcomplete)/100) %>%
  group_by(Birth_year, Date_of_statistics) %>%
  summarise(nrvacc = sum(Populatie * covprimcomplete),
            Populatie = sum(Populatie)) %>%
  group_by(Date_of_statistics) %>%
  summarise(propvacc = max(nrvacc) / popsize) %>%
  ungroup() %>%
  full_join(
    tibble(Date_of_statistics = seq(as.Date("2020-02-01"), max(vaccdata$Date_of_statistics), 1))
  ) %>%
  arrange(Date_of_statistics) %>%
  mutate(propvacc = if_else(Date_of_statistics == min(Date_of_statistics), 0, propvacc),
         propvacc = zoo::na.approx(propvacc),
         propprot = propvacc * 0.75) %>%
  rename(date = Date_of_statistics)

###################
### Google mobility
###################
## source: https://www.google.com/covid19/mobility/
# # Google mobility report
googlemobility <- bind_rows(
  read_csv(file = "data/2020_NL_Region_Mobility_Report.csv") %>%
    filter(is.na(sub_region_1)),
  read_csv(file = "data/2021_NL_Region_Mobility_Report.csv") %>%
    filter(is.na(sub_region_1)),
  read_csv(file = "data/2022_NL_Region_Mobility_Report.csv") %>%
    filter(is.na(sub_region_1)))

googlemobility <- googlemobility %>%
  rename(work = workplaces_percent_change_from_baseline,
         grocery = grocery_and_pharmacy_percent_change_from_baseline,
         residential = residential_percent_change_from_baseline,
         transit = transit_stations_percent_change_from_baseline,
         retail = retail_and_recreation_percent_change_from_baseline,
         parks = parks_percent_change_from_baseline) %>%
  dplyr::select(date, retail:residential)

########################
### Intervention periods
########################
# Intervention periods
period_dates <- as.Date(c("2020-01-01","2020-03-13", "2020-03-28", "2020-05-10", "2020-06-01",
                          "2020-07-05", "2020-08-30", "2020-09-28", "2020-10-14",
                          "2020-10-25", "2020-11-04", "2020-11-18", "2020-12-14",
                          "2020-12-20", "2021-01-03", "2021-01-22", "2021-02-07",
                          "2021-02-28", "2021-04-18", "2021-04-27", "2021-05-02",
                          "2021-05-16", "2021-05-18", "2021-05-30", "2021-06-04",
                          "2021-06-25", "2021-07-09", "2021-08-08", 
                          "2021-08-31", "2021-09-24", "2021-10-17",
                          "2021-10-31", "2021-11-12", "2021-11-27",
                          "2021-12-18", "2021-12-24", "2022-01-09",
                          "2022-01-16", "2022-01-25", "2022-02-17", "2022-02-24", "2022-04-01"))
period_labels <- paste0("_", substr(seq(101, 99 + length(period_dates)), 2, 3))

intervention_periods <- tibble(date = seq(as.Date("2020-01-01"), as.Date("2022-04-01"), by = "day")) %>% 
  mutate(period = cut(date,
                      breaks = period_dates + 1,
                      labels = period_labels))
rm(period_dates, period_labels)

###############
### KNMI-data
#############
dataKNMI_raw <- read_csv(file = "data/KNMI_20230711.txt", col_names = TRUE, comment = "#") %>%
  dplyr::select(STN, YYYYMMDD, TG, TN, TX, UG)
names(dataKNMI_raw) <- c("station", "date", "temp_avg", "temp_min", "temp_max", "rh_avg")

# from Wu et al (2012), https://academic.oup.com/jid/article/206/12/1862/865907
absolute_humidity <- function(rel_hum, tempC) rel_hum*6.112*2.1674*exp(17.67*tempC/(243.5+tempC))/(273.15+tempC)

dataKNMI <- dataKNMI_raw %>% mutate(date = as.Date(as.character(date), format = "%Y%m%d"),
                                    ah = absolute_humidity(rel_hum = rh_avg, tempC = temp_avg/10),
                                    temp = temp_avg/10,
                                    dayinyear = format(date, "%d-%m")) %>% 
  filter(date >= as.Date("2000-01-01")) %>%
  mutate(include4mean = date <= ymd("2019-12-31"),
         ah_reference = mean(ah[include4mean]),
         temp_reference = mean(temp[include4mean]))  %>%
  mutate(ah = ah - ah_reference,
         temp = temp - temp_reference) %>%
  group_by(dayinyear) %>% 
  mutate(ah_avg = mean(ah[include4mean]),
         temp_avg = mean(temp[include4mean])) %>%
  ungroup() %>%
  dplyr::select(-include4mean,-temp_min,-temp_max,-rh_avg,-station) %>%
  mutate(date_04 = date + 4*365 + 1,
         date_08 = date + 8*365 + 2,
         date_12 = date + 12*365 + 3,
         date_16 = date + 16*365 + 4,
         date_20 = date + 20*365 + 5) %>%
  rename(date_00 = date) %>%
  pivot_longer(c(date_00, date_04:date_20),
               names_to = "yearsago",
               values_to = "date") %>%
  filter(date < max(date) - 20*365 - 4) %>%
  mutate(yearsago = substring(yearsago, 6)) %>%
  pivot_wider(names_from = yearsago, values_from = ah:temp) %>%
  rename(ah = ah_00,
         temp = temp_00) %>%
  filter(date >= as.Date("2020-01-01"))
rm(dataKNMI_raw)

#####################################
### Nederlands Verplaatsingspanel ###
### https://www.goudappel.nl/nl/expertises/data-en-it-oplossingen/nederlands-verplaatsingspanel/het-nvp-tijdens-covid-19
### Accessed: 2023-09-15
#####################################
NVPdata <- read_csv("data/NVP.csv")
NVPdata <- NVPdata %>%
  slice(rep(1:n(), each = 7)) %>%
  mutate(date = seq.Date(from = as.Date("2020-03-02"), by = 1, length.out = 7*182)) %>%
  dplyr::select(-concat, -maand)

#################################
### combine all data and add lags
#################################
alldata_daily <- epicurve %>% 
  mutate(logRt = lead(log(instR), 5),
         weight = lead(1/instRlogvar, 5)) %>%
  dplyr::select(date, logRt, weight, dayinweek) %>%
  full_join(variants) %>%
  full_join(vaccdata) %>%
  full_join(susceptibles) %>%
  mutate(susceptibles = susceptibles/popsize,
         susceptibles = if_else(is.na(susceptibles), min(susceptibles, na.rm = TRUE), susceptibles),
         unprotecteds = (1 - propprot) * (susceptibles))  %>%
  full_join(dataKNMI) %>% 
  full_join(googlemobility) %>%
  arrange(date) %>%
  mutate(across(retail:residential, ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE))) %>%
  full_join(intervention_periods) %>%
  filter(date >= as.Date("2020-02-15") & date <= as.Date("2022-04-01")) %>%
  mutate(logRvariant = log(variantcurve),
         logUnprotected = log(unprotecteds)) %>%
  filter(!is.na(logRt)) %>%
  dplyr::select(-(variantcurve:unprotecteds))

# weekly averages of weather and mobility variables in week data
alldata_weekly <- alldata_daily %>%
  mutate(across(ah:residential, ~ zoo::rollmean(.x, k = 7, fill = "extend"))) %>%
  left_join(NVPdata, by = "date")
alldata_weekly <- alldata_weekly %>%
  mutate(across(starts_with("index"), ~if_else(is.na(.x), 100, .x))) %>%
  mutate(across(starts_with("index"), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE))) 
  
# make lagged variables in daily data
for(i in 1:6) alldata_daily <- alldata_daily %>% mutate(across(retail:residential, ~lag(.x, i), .names = paste0("{.col}_",i)))
for(i in 1:7) alldata_daily <- alldata_daily %>% mutate(across(c(ah,temp), ~lag(.x, i), .names = paste0("{.col}_",i)))
for(i in 1:7) alldata_daily <- alldata_daily %>% mutate(across(logRt, ~lag(.x, i), .names = paste0("{.col}_",i)))

# make lagged variables in weekly data
alldata_weekly <- alldata_weekly %>% 
  group_by(dayinweek) %>% 
  mutate(across(c(ah, temp, ah_avg, temp_avg, retail:residential, index_thuisblijvers:index_zorg), ~lag(.x, 1), .names = paste0("{.col}_",1))) %>% 
  ungroup()
alldata_weekly <- alldata_weekly %>% 
  group_by(dayinweek) %>% 
  mutate(across(logRt, ~lag(.x, 1), .names = paste0("{.col}_",1))) %>% 
  ungroup()

# remove data from early February, which were used to make the lagged variables
alldata_daily <- alldata_daily %>%
  filter(!is.na(logRt_7))
alldata_weekly <- alldata_weekly %>%
  filter(!is.na(logRt_1))

# split "dansen met Janssen" period into two, and create new datasets
alldata_daily_splitJanssen <- alldata_daily %>% 
  mutate(period = if_else(period == "_26" & logRt > 0.5, "_26J", period))
alldata_weekly_splitJanssen <- alldata_weekly %>% 
  mutate(period = if_else(period == "_26" & logRt > 0.5, "_26J", period))

##################################################
### Oxford Stringency Index, only for plotting ###
##################################################
OxfordIndex <- read_csv("data/OxCGRT_timeseries_StringencyIndex_v1.csv") %>%
  filter(CountryName == "Netherlands") %>%
  dplyr::select(-(1:7)) %>%
  pivot_longer(everything(), names_to = "date", values_to = "OSI") %>%
  mutate(date = as.Date(date, "%d-%b-%y")) %>%
  filter(date > ymd("2020-02-22") & date < ymd("2022-04-01"))


rm(dataKNMI, epicurve, googlemobility, intervention_periods, NVPdata, 
   susceptibles, vaccdata, variants, i, popsize, absolute_humidity)

