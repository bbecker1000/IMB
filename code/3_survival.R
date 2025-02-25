library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)
library(ggplot2)
library(survival)
library(icenReg)

data <- read_csv(here::here("data", "dd_data.csv")) %>% 
  mutate(
    imb_id = as.factor(imb_id),
    host_plant = as.factor(host_plant),
    overall_survival = as.factor(overall_survival),
    stage = as.factor(stage),
    start = ymd(start),
    end = ymd(end)
  ) %>% 
  mutate_at(vars(contains("stage")), 
            factor,
            levels = c("egg", "first", "second", "third", "fourth", "fifth", "pupa"),
            ordered = TRUE) %>% 
  mutate(
    event_type = if_else(survival == TRUE, 0, 3)
  ) # need to make times numeric....

stage1 <- data %>% 
  filter(stage == "first")

stage2 <- data %>% 
  filter(stage == "second")

stage3 <- data %>% 
  filter(stage == "third")


# here, the "event" is death, and "lack of event" is progression to the next stage
# we know the exact time of progression to the next stage, but only an interval in which death occurred
# for this, we would put each stage into a separate model
# and look at survival separately in each stage

m.stage1 <- ic_sp(
  Surv(time = start, time2 = end, event = event_type, type = 'interval2') ~ 1,
  data = stage1
)
