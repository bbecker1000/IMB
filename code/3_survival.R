library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)
library(ggplot2)
library(survival)
library(icenReg)
library(flexsurv)

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
            levels = c("first", "second", "third", "fourth", "fifth", "pupa"),
            ordered = TRUE) %>% 
  mutate(
    event = as.factor(if_else(survival, paste0("progression", (as.integer(stage) + 1)), "death")),
    status = 1
  ) %>% 
  filter(duration > 0)

write_csv(data, "data.csv")

# need to add rows representing start times for adult and death
# data2 <- data %>% 
#   group_by(imb_id) %>%
#   slice_max(stage) %>% 
#   mutate( # 8 = death, 7 = adult, 9 = censored
#     state = if_else(survival == FALSE, 8, if_else(state == 6, 7, 9)),
#     start = end,
#     end = NA
#   )

# trying to use flexsurv, hopefully will be able to just to one model

fs.st1 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                    data = data %>% filter(stage == "first") %>% droplevels(.),
                    dists = c(progression2 = "weibull", death = "exponential"))

fs.st2 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data %>% filter(stage == "second") %>% droplevels(.),
                      dists = c(progression3 = "weibull", death = "exponential"))

fs.st3 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data %>% filter(stage == "third") %>% droplevels(.),
                      dists = c(progression4 = "weibull", death = "exponential"))

fs.st4 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data %>% filter(stage == "fourth") %>% droplevels(.),
                      dists = c(progression5 = "weibull", death = "exponential"))

fs.st5 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data %>% filter(stage == "fifth") %>% droplevels(.),
                      dists = c(progression6 = "weibull", death = "exponential"))

fs.stp <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data %>% filter(stage == "pupa") %>% droplevels(.),
                      dists = c(progression7 = "weibull", death = "exponential"))

fs.total <- fmixmsm("1st instar" = fs.st1, "2nd instar" = fs.st2, "3rd instar" = fs.st3)

ppath_fmixmsm(fs.total, B = 100)

qfinal_fmixmsm(fs.total, B = 10)

ajfit_flexsurvmix(fs.st2, startname = "second")

# trying to use icenReg, requires separate models for each state
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
