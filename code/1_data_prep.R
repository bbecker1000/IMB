library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)
library(climate)

setwd(here::here("code"))

# reading in rearing data, assigning data types
raw_data <- read_csv(here::here("data", "IMB_RearingReleaseData.csv")) %>% 
  rename(imb_id = `IMB ID`,
         rearing_year = `REARING YEAR`,
         host_plant = `HOST PLANT`,
         collection_site = `COLLECTION SITE`,
         collection_area = `COLLECTION AREA`,
         collection_date = `COLLECTION DATE`,
         collection_stage = `COLLECTION STAGE`,
         hatch_date = `HATCH DATE`,
         date_instar_2 = `Date of 2nd instar`,
         date_instar_3 = `Date 3rd instar`,
         date_instar_4 = `Date 4th instar`,
         date_instar_5 = `Date fifth instar`,
         pupation_date = `PUPATION DATE`,
         larval_days = `LARVAL DAYS`,
         overall_survival = `SURVIVE TO PUPATION`,
         stage_at_death = `STAGE AT DEATH`,
         notes = `REARING NOTES`,
         release_year = `RELEASE YEAR`,
         eclosure_date = `ECLOSURE DATE`,
         release_date = `RELEASE DATE`,
         release_site = `RELEASE SITE`,
         release_area = `RELEASE AREA`,
         sex = SEX,
         lat = `LAT.`,
         long = `LONG.`,
         release_notes = `RELEASE NOTES`
         ) %>% 
  mutate(larval_days = if_else((larval_days == "#VALUE!" | larval_days == "#REF!"), NA, larval_days),
         overall_survival = if_else(overall_survival == "y", "Y", overall_survival),
         overall_survival = if_else((overall_survival == "n" | overall_survival == "NV"), "N", overall_survival),
         sex = if_else(sex == "f", "F", sex),
         sex = if_else(sex == "m", "M", sex),
         sex = if_else(sex == "UNK", NA, sex),
         rearing_year = as.integer(rearing_year),
         host_plant = as.factor(host_plant),
         collection_site = as.factor(collection_site),
         collection_area = as.factor(collection_area),
         collection_stage = as.factor(collection_stage),
         larval_days = as.integer(larval_days),
         overall_survival = as.factor(overall_survival),
         stage_at_death = as.factor(stage_at_death),
         release_year = as.integer(release_year),
         release_site = as.factor(release_site),
         release_area = as.factor(release_area),
         sex = as.factor(sex)) %>% 
  mutate_at(vars(contains("date")), mdy) %>% 
  mutate(collection_stage = fct_collapse(collection_stage, 
                                       egg = c("egg", "Egg"),
                                       first = c("First", "1st"),
                                       second = c("Second", "2nd"),
                                       third = c("Third", "3rd"),
                                       fourth = c("Fourth", "4th"),
                                       fifth = c("Fifth"),
                                       other_level = NA
                                       ),
         stage_at_death = fct_collapse(stage_at_death,
                                       egg = c("Egg"),
                                       first = c("First", "1st"),
                                       second = c("Second", "2nd"),
                                       third = c("Third", "3rd"),
                                       fourth = c("Fourth"),
                                       fifth = c("Fifth", "5th"),
                                       pupa = c("pupa", "Pupa", "pupae")
         ))

# reading in temperature data
threshold <- 10 # really not sure what the temperature threshold is supposed to be...
# cannot find minimum development threshold for IMB, but i found a source that generally says the min threshold for lots
# of butterflies is 15

# adding temperature data from weather station

station_temps <- meteo_noaa_hourly(
  station = "727985-94276",
  year = 2013:2024,
  fm12 = FALSE
) 

station_temps_cleaned <- station_temps %>% 
  select(date, year, month, day, hour, t2m) %>% 
  rename(temp = t2m) %>% 
  mutate(
    date = ymd(substr(date, 1, 10))
  ) %>% 
  group_by(date) %>% 
  summarise(
     mean_temp = mean(temp, na.rm = TRUE)
     # max_temp = max(temp, na.rm = TRUE),
     # min_temp = min(temp, na.rm = TRUE),
     # avg_temp = (max_temp + min_temp)/2
  ) %>% 
  ungroup() %>% 
  mutate(degree_days = pmax(0, mean_temp - threshold)) %>% 
  select(-mean_temp)


# removing illogical data, adding durations for analysis, rearranging data to make analysis easier
cleaned_data <- raw_data %>% 
  select(-notes, -release_notes, -release_site, -release_area, -lat, -long) %>% 
  mutate(larval_days = if_else(larval_days < 0 | larval_days > 100, NA, larval_days)) %>% 
  filter(!is.na(collection_date),
         is.na(hatch_date) | collection_date <= hatch_date,
         is.na(hatch_date) | is.na(date_instar_2) | hatch_date <= date_instar_2,
         is.na(date_instar_2) | is.na(date_instar_3) | date_instar_2 <= date_instar_3,
         is.na(date_instar_3) | is.na(date_instar_4) | date_instar_3 <= date_instar_4,
         is.na(date_instar_4) | is.na(date_instar_5) | date_instar_4 <= date_instar_5,
         is.na(date_instar_5) | is.na(pupation_date) | date_instar_5 <= pupation_date,
         !(is.na(hatch_date) & is.na(date_instar_2) & is.na(date_instar_3) & is.na(date_instar_4) & is.na(date_instar_5) & is.na(pupation_date))) %>% 
  mutate(
    duration_first = interval(hatch_date, date_instar_2),
    duration_second = interval(date_instar_2, date_instar_3),
    duration_third = interval(date_instar_3, date_instar_4),
    duration_fourth = interval(date_instar_4, date_instar_5),
    duration_fifth = interval(date_instar_5, pupation_date),
    duration_pupa = interval(pupation_date, eclosure_date)
  ) %>% 
  mutate(
    overall_survival = if_else(is.na(eclosure_date), "N", overall_survival),
    stage_at_death = if_else((is.na(eclosure_date) & is.na(stage_at_death)), "pupa", stage_at_death)
  )

# rearranging data to make analysis easier
data <- cleaned_data %>% 
  select(-contains("date"), -release_year, -collection_site, -collection_area, -sex) %>% 
  pivot_longer(cols = starts_with("duration"), names_to = "stage", names_prefix = "duration_", values_to = "duration") %>% 
  mutate_at(vars(contains("stage")), 
            factor,
            levels = c("egg", "first", "second", "third", "fourth", "fifth", "pupa"),
            ordered = TRUE) %>% 
  mutate(invalid_stage = (collection_stage > stage) | stage_at_death <= stage,
         invalid_stage = if_else(is.na(invalid_stage), FALSE, invalid_stage)) %>% 
  filter(invalid_stage == FALSE,
         !is.na(duration),
         !is.na(overall_survival)) %>%  
  mutate(survival = if_else(as.character(overall_survival) == "Y", TRUE, (as.integer(stage) + 1) != as.integer(stage_at_death))) %>% 
  select(-invalid_stage, -stage_at_death, -larval_days) %>% 
  mutate(
    start = int_start(duration),
    end = int_end(duration),
  )


data_no_temp <- data %>% 
  mutate(
    duration = as.integer((end-start) / 86400) #value in days
  )

write_csv(data_no_temp, here::here("data", "data_no_temp.csv"))

# adding degree day data to rearing data
result <- station_temps_cleaned %>%
  cross_join(data) %>%  # Cross join
  filter(date >= start & date <= end) %>% # Keep only valid date ranges
  group_by(imb_id, stage, start, end) %>%        # Group by unique intervals
  summarise(total_degree_days = sum(degree_days), .groups = "drop")

# Merge summed results back into the original data
dd_data <- data %>%
  left_join(result, by = c("imb_id", "stage", "start", "end")) %>% 
  filter(!is.na(total_degree_days)) %>% 
  mutate(
    duration = as.integer((end-start) / 86400) #value in days
  )


write_csv(dd_data, here::here("data", "dd_data.csv"))
