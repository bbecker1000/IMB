library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)

setwd(here::here("code"))

# reading in data, assigning data types
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
         survival = `SURVIVE TO PUPATION`,
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
         survival = if_else(survival == "y", "Y", survival),
         survival = if_else((survival == "n" | survival == "NV"), "N", survival),
         sex = if_else(sex == "f", "F", sex),
         sex = if_else(sex == "m", "M", sex),
         sex = if_else(sex == "UNK", NA, sex),
         rearing_year = as.integer(rearing_year),
         host_plant = as.factor(host_plant),
         collection_site = as.factor(collection_site),
         collection_area = as.factor(collection_area),
         collection_stage = as.factor(collection_stage),
         larval_days = as.integer(larval_days),
         survival = as.factor(survival),
         stage_at_death = as.factor(stage_at_death),
         release_year = as.integer(release_year),
         release_site = as.factor(release_site),
         release_area = as.factor(release_area),
         sex = as.factor(sex)) %>% 
  mutate_at(vars(contains("date")), mdy)

# removing illogical data
cleaned_data <- raw_data %>% 
  mutate(larval_days = if_else(larval_days < 0 | larval_days > 100, NA, larval_days)) %>% 
  filter(!is.na(collection_date),
         is.na(hatch_date) | collection_date <= hatch_date,
         is.na(hatch_date) | hatch_date <= date_instar_2,
         is.na(date_instar_2) | date_instar_2 <= date_instar_3,
         is.na(date_instar_3) | date_instar_3 <= date_instar_4,
         is.na(date_instar_4) | date_instar_4 <= date_instar_5,
         is.na(date_instar_5) | date_instar_5 <= pupation_date)

# getting a sense of NA values in each column

raw_data %>% level(survival)
raw_data %>% count(sex)
raw_data %>% count(larval_days)

raw_data %>% count(complete.cases(.))

diff <- anti_join(raw_data, cleaned_data, by = join_by(imb_id))

weird_data <- raw_data %>% 
  filter((!is.na(hatch_date) & collection_date > hatch_date) |
         (!is.na(hatch_date) & hatch_date > date_instar_2) |
         (!is.na(date_instar_2) & date_instar_2 > date_instar_3) |
         (!is.na(date_instar_3) & date_instar_3 > date_instar_4) |
         (!is.na(date_instar_4) & date_instar_4 > date_instar_5) |
         (!is.na(date_instar_5) & date_instar_5 > pupation_date))


