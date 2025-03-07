library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)
library(climate)

setwd(here::here("code"))

#### reading in rearing data, assigning data types ####
raw_data <- read_csv(here::here("data", "IMB_RearingReleaseData.csv")) %>% 
  rename(imb_id = `IMB ID`, rearing_year = `REARING YEAR`, host_plant = `HOST PLANT`, collection_site = `COLLECTION SITE`,
  collection_area = `COLLECTION AREA`, collection_date = `COLLECTION DATE`, collection_stage = `COLLECTION STAGE`,
  hatch_date = `HATCH DATE`, date_instar_2 = `Date of 2nd instar`, date_instar_3 = `Date 3rd instar`,
  date_instar_4 = `Date 4th instar`, date_instar_5 = `Date fifth instar`, pupation_date = `PUPATION DATE`,
  larval_days = `LARVAL DAYS`, overall_survival = `SURVIVE TO PUPATION`, stage_at_death = `STAGE AT DEATH`,
  notes = `REARING NOTES`, release_year = `RELEASE YEAR`, eclosure_date = `ECLOSURE DATE`,
  release_date = `RELEASE DATE`, release_site = `RELEASE SITE`, release_area = `RELEASE AREA`,
  sex = SEX, lat = `LAT.`, long = `LONG.`, release_notes = `RELEASE NOTES`) %>% 
  select(-notes, -release_notes, -release_site, -release_area, -lat, -long, -collection_site, -collection_area, -larval_days) %>% 
  mutate(overall_survival = if_else(overall_survival == "y", "Y", overall_survival),
         overall_survival = if_else((overall_survival == "n" | overall_survival == "NV"), "N", overall_survival),
         sex = if_else(sex == "f", "F", if_else(sex == "m", "M", if_else(sex == "UNK", NA, sex))),
         host_plant = as.factor(host_plant),
         collection_stage = as.factor(collection_stage),
         overall_survival = as.factor(overall_survival),
         stage_at_death = as.factor(stage_at_death),
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

##### removing illogical data, adding durations for analysis, rearranging data to make analysis easier ####
cleaned_data <- raw_data %>% 
  filter(collection_stage == "egg",
        !is.na(collection_date),
         is.na(hatch_date) | collection_date <= hatch_date,
         is.na(hatch_date) | is.na(date_instar_2) | hatch_date <= date_instar_2,
         is.na(date_instar_2) | is.na(date_instar_3) | date_instar_2 <= date_instar_3,
         is.na(date_instar_3) | is.na(date_instar_4) | date_instar_3 <= date_instar_4,
         is.na(date_instar_4) | is.na(date_instar_5) | date_instar_4 <= date_instar_5,
         is.na(date_instar_5) | is.na(pupation_date) | date_instar_5 <= pupation_date,
         is.na(pupation_date) | is.na(eclosure_date) | pupation_date <= eclosure_date,
         !(is.na(date_instar_2) & is.na(date_instar_3) & is.na(date_instar_4) & is.na(date_instar_5))) %>% 
  mutate(
    # duration_egg = interval(collection_date, hatch_date),
    duration_first = interval(hatch_date, date_instar_2),
    duration_second = interval(date_instar_2, date_instar_3),
    duration_third = interval(date_instar_3, date_instar_4),
    duration_fourth = interval(date_instar_4, date_instar_5),
    duration_fifth = interval(date_instar_5, pupation_date),
    duration_pupa = interval(pupation_date, eclosure_date),
    overall_survival = if_else(is.na(eclosure_date), "N", overall_survival),
    final_stage = case_when(
      !is.na(eclosure_date) ~ "adult",
      !is.na(pupation_date) ~ "pupa",
      !is.na(date_instar_5) ~ "fifth",
      !is.na(date_instar_4) ~ "fourth",
      !is.na(date_instar_3) ~ "third",
      !is.na(date_instar_2) ~ "second",
      !is.na(hatch_date) ~ "first",
      TRUE ~ "egg"
      )
    ) 

data <- cleaned_data %>% 
  select(-contains("date"), -release_year, -stage_at_death) %>% 
  pivot_longer(
    cols = starts_with("duration"), 
    names_to = "stage", 
    names_prefix = "duration_", 
    values_to = "duration"
    ) %>% 
  mutate(
    stage = factor(stage, levels = c("egg", "first", "second", "third", "fourth", "fifth", "pupa", "adult"), ordered = TRUE),
    start = int_start(duration),
    end = int_end(duration),
    survival = !is.na(end)
  ) %>% 
  filter(!is.na(start)) %>%
  arrange(imb_id, stage) %>% 
  group_by(imb_id) %>% 
  mutate(
    missing_stage = is.na(end),
    later_all_na = cumall(rev(cumall(rev(is.na(end))))),
    status = case_when(
      !missing_stage ~ "complete",
      (missing_stage & later_all_na) | (stage == final_stage) ~ "died",
      missing_stage & !later_all_na ~ "censored"
    )
  )

data_no_censor <- data %>%
  group_by(imb_id) %>% 
  filter(!("censored" %in% status) | status == "complete") %>% 
  ungroup()

# censored_individuals <- data %>% 
#   filter(status == "censored") %>% 
#   group_by(imb_id) %>% 
#   arrange(stage) %>% 
#   slice_head()

data <- data_no_censor %>% 
  arrange(imb_id) %>% 
  select(-overall_survival, -missing_stage, -later_all_na) %>% 
  group_by(imb_id) %>% 
  mutate(d_dups = sum(status == "died")) %>% 
  filter(d_dups <= 1) %>% 
  select(-d_dups, -survival)


mean_durations <- data %>% 
  filter(!is.na(end)) %>% 
  group_by(stage) %>% 
  summarise(avg_length = mean(end-start))

data <- left_join(data, mean_durations, by = c("stage")) %>% 
  mutate(end = if_else(status != "died", end, round_date(start + avg_length, unit = "day")),
         stage = if_else(status == "died", "death", stage),
         stage = factor(stage, levels = c("egg", "first", "second", "third", "fourth", "fifth", "pupa", "adult", "death"), ordered = TRUE),
         status = as_factor(status),
         start = if_else(start == end, start - ddays(1), start)) %>% 
  select(-avg_length) %>% 
  group_by(imb_id) %>% 
  mutate(collection_day = min(start)) %>% 
  ungroup()


#### temperature data ####  
# adding temperature data from weather station

# reading in temperature data
threshold <- 10 # really not sure what the temperature threshold is supposed to be...
# cannot find minimum development threshold for IMB, but i found a source that generally says the min threshold for lots
# of butterflies is 15

# only run this once, it takes a while. After that just read_csv to get the local data
# station_temps <- meteo_noaa_hourly(
#   station = "727985-94276",
#   year = 2013:2024,
#   fm12 = FALSE
# ) 
# write_csv(station_temps, here::here("data", "station_temps.csv"))

station_temps <- read_csv(here::here("data", "station_temps.csv"))

station_temps_cleaned <- station_temps %>% 
  select(date, year, month, day, hour, t2m) %>% 
  rename(temp = t2m) %>% 
  mutate(
    date = ymd(substr(date, 1, 10))
  ) %>% 
  group_by(date) %>% 
  summarise(
    mean_temp = mean(temp, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(degree_days = pmax(0, mean_temp - threshold)) %>% 
  select(-mean_temp)

#### combining data with temp data ####
# adding degree day data to rearing data
result <- station_temps_cleaned %>%
  cross_join(data) %>%  # Cross join
  filter(date >= collection_day & date <= end) %>% # Keep only valid date ranges
  group_by(imb_id, stage, start, end) %>%        # Group by unique intervals
  summarize(cum_degree_days = sum(degree_days, na.rm = TRUE)) %>% 
  ungroup()
  # group_by(imb_id, stage, collection_day, end) %>% 
  # mutate(total_degree_days = sum(degree_days, na.rm = TRUE))

# Merge summed results back into the original data
dd_data <- data %>%
  left_join(result, by = c("imb_id", "stage", "start", "end")) %>% 
  mutate(
    duration = as.integer((end-start))
  ) %>% 
  mutate(
    # this is where I should do starting and ending dd too
    start = as.integer(start - collection_day) / 86400, # seconds to days
    end = as.integer(start + duration),
    start = if_else(end == start, start - 1, start),
    stage = factor(stage, levels = c("first", "second", "third", "fourth", "fifth", "pupa", "adult", "death"), ordered = TRUE),
    stage2 = fct_recode(stage, "second" = "first", "third" = "second", "fourth" = "third", "fifth" = "fourth", "pupa" = "fifth", "adult" = "pupa"),
    after_2018 = if_else(rearing_year >= 2018, TRUE, FALSE)
  ) %>% 
  select(-final_stage, -status, -collection_day) %>% 
  group_by(imb_id) %>% 
  arrange(imb_id, stage) %>% 
  mutate(
    start_dd = if_else(row_number() == 1, 0, lag(cum_degree_days))
    # total_dd = start_dd + degree_days
  ) %>% 
  ungroup() %>% 
  filter(imb_id != "23SI024",
         imb_id != "18SI029")


dd_data <- dd_data %>% 
  group_by(imb_id) %>% 
  mutate(to_compare = row_number() == 1,
  to_remove = if_else(to_compare & stage != "first", "yes", "no"),
  to_remove = if_else("yes" %in% to_remove, T, F)) %>% 
  filter(to_remove == FALSE) %>% 
  select(-to_compare, -to_remove)
  
write_csv(dd_data, here::here("data", "dd_data.csv"))
