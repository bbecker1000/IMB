library(dplyr)
library(tidyverse)
library(here)

setwd(here::here("Code"))

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
         sex = if_else(sex == "UNK", NA, sex))

# getting a sense of NA values in each column

raw_data %>% count(survival)
raw_data %>% count(sex)
raw_data %>% count(larval_days)

