library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)
library(ggplot2)
library(msm)

# this is not going well right now :(

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
            ordered = TRUE)

# need to convert dates to days since start of observations

start_date <- min(data$start, na.rm = TRUE)

msm_data <- data %>%
  mutate(state = as.numeric(factor(stage, levels = c("first", "second", "third", "fourth", "fifth", "pupa")))) %>%
  select(imb_id, start, end, state, survival) %>%
  arrange(imb_id, state) %>% 
  mutate(
    start = as.numeric(difftime(start, start_date, units = "days")),
    end = as.numeric(difftime(end, start_date, units = "days"))
  )

# need to add rows representing start times for adult and death
msm_data2 <- msm_data %>% 
  group_by(imb_id) %>%
  slice_max(state) %>% 
  mutate( # 8 = death, 7 = adult, 9 = censored
    state = if_else(survival == FALSE, 8, if_else(state == 6, 7, 9)),
    start = end,
    end = NA
  )

msm_data_all_states <- full_join(msm_data, msm_data2) %>%
  group_by(imb_id, start) %>% 
  mutate(n = n()) %>% 
  filter(n == 1) %>% 
  select(imb_id, start, state) %>% 
  ungroup() %>% 
  arrange(imb_id, start)


statetable.msm(state, imb_id, data = msm_data_all_states)

Q <- rbind(
  # 1  2   3   4   5   p   a  d
  c( 0, .75,  0,  0,  0,  0,  0, .1, .1),  # First → Second OR First → Death OR → Censored
  c( 0, 0,  .75,  0,  0,  0,  0, .1, .1),  # Second → Third OR Second → Death OR → Censored
  c( 0, 0,  0,  .75,  0,  0,  0, .1, .1),  # Third → Fourth OR Third → Death OR → Censored
  c( 0, 0,  0,  0,  .75,  0,  0, .1, .1),  # Fourth → Fifth OR Fourth → Death OR → Censored
  c( 0, 0,  0,  0,  0,  .75,  0, .1, .1),  # Fifth → Pupa OR Fifth → Death OR → Censored
  c( 0, 0,  0,  0,  0,  0,  .75, .1,  0),  # Pupa → Adult OR Pupa → Death OR → Censored
  c( 0, 0,  0,  0,  0,  0,  0, 0, 0),  # Adult (absorbing state)
  c( 0, 0,  0,  0,  0,  0,  0, 0, 0),  # Death (absorbing state)
  c( 0, 0,  0,  0,  0,  0,  0, 0, 0)   # Censored (absorbing state)
)

Q.simple <- rbind(
  # 1     2   3   4   5   p   a    d
  c( 0, .75,  0,  0,  0,  0,  0, .1),  # First → Second OR First → Death
  c( 0, 0,  .75,  0,  0,  0,  0, .1),  # Second → Third OR Second → Death
  c( 0, 0,  0,  .75,  0,  0,  0, .1),  # Third → Fourth OR Third → Death 
  c( 0, 0,  0,  0,  .75,  0,  0, .1),  # Fourth → Fifth OR Fourth → Death
  c( 0, 0,  0,  0,  0,  .75,  0, .1),  # Fifth → Pupa OR Fifth → Death
  c( 0, 0,  0,  0,  0,  0,  .75, .1),  # Pupa → Adult OR Pupa → Death
  c( 0, 0,  0,  0,  0,  0,  0, 0),  # Adult (absorbing state)
  c( 0, 0,  0,  0,  0,  0,  0, 0)  # Death (absorbing state)
)

Q.crude <- crudeinits.msm(state ~ start, imb_id, data = msm_data_all_states, censor = 9, qmatrix = Q.simple)


msm_fit <- msm(state ~ start, subject = imb_id, data = msm_data_all_states,
               qmatrix = Q.crude, gen.inits = F, obstype = 1, opt.method = "optim", censor = 9,
               control = list(maxit = 5000))
msm_fit
