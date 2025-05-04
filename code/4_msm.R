library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)
library(ggplot2)
library(msm)

# keeping this file as an archive, not using any multi-state models for this right now

#### trying to use flexsurv, hopefully will be able to just to one model ####

fs.st1 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data_fs_2 %>% filter(stage == "first"),
                      dists = c(progression = "weibull", death = "exponential"))

fs.st2 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data_fs_2 %>% filter(stage == "second"),
                      dists = c(progression = "weibull", death = "exponential"))

fs.st3 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data_fs_2 %>% filter(stage == "third"),
                      dists = c(progression = "weibull", death = "exponential"))

fs.st4 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data_fs_2 %>% filter(stage == "fourth"),
                      dists = c(progression = "weibull", death = "exponential"),
                      fixedpars = TRUE,
                      method = "direct")

fs.st5 <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data_fs_2 %>% filter(stage == "fifth"),
                      dists = c(progression = "weibull", death = "exponential"))

fs.stp <- flexsurvmix(Surv(duration, status) ~ 1, event = event,
                      data = data_fs_2 %>% filter(stage == "pupa"),
                      dists = c(progression = "weibull", death = "exponential"))

fs.total <- fmixmsm("1st instar" = fs.st1, "2nd instar" = fs.st2,  "3rd instar" = fs.st3)

ppath_fmixmsm(fs.total, B = 100, final = FALSE)

qfinal_fmixmsm(fs.total, B = 10)

p_flexsurvmix(fs.st1, startname = "first instar", newdata = c(seq(0, 15, by = 1)))

plot_data <- ajfit_flexsurvmix(fs.st1, startname = "first instar") %>% 
  filter(model == "Parametric mixture")

plot_data <- ajfit_flexsurvmix(fs.st2, startname = "second instar") %>% 
  filter(model == "Parametric mixture")

plot_data <- ajfit_flexsurvmix(fs.st3, startname = "third instar") %>% 
  filter(model == "Parametric mixture")

ggplot(plot_data, aes(x = time, y = val)) +
  geom_line(aes(color = state)) +
  theme_bw()




#### flexsurv framework 1 ####

data_fs_1 <- data %>% 
  select(-rearing_year, -host_plant, -overall_survival) %>% 
  mutate(
    stage = as.integer(stage),
    duration = if_else(duration == 0, 0.5, duration)
  ) %>% 
  filter(stage == 1 | stage == 2 | stage == 3) %>% 
  rename(from = stage) %>% 
  mutate(
    to = if_else(survival, from + 1, 0)
  ) %>% 
  select( -survival) %>%
  group_by(imb_id) %>% 
  mutate(
    collection_day = min(start)
  ) %>% 
  ungroup() %>% 
  mutate(
    start = as.integer(start - collection_day),
    end = as.integer(start + duration),
    status = 1, # observed
    trans = paste0(from, to)
  ) %>% 
  relocate(
    imb_id, from, to, start, end, status, trans, duration, collection_day
  )

data_fs_1_censored <- data_fs_1 %>% 
  mutate(to = if_else(to == 0, from + 1, 0),
         status = 0, #censored
         trans = paste0(from, to)
  ) 

data_fs_1 <- rbind(data_fs_1, data_fs_1_censored) %>% 
  mutate(trans = as.factor(trans))

# generates curvy curves
m2.fs1 <- flexsurvreg(Surv(duration, status) ~ trans + shape(trans), data = data_fs_1, 
                      dist = "weibull")

# doesn't work
m3.fs1 <- flexsurvreg(Surv(start, end, status) ~ trans + shape(trans), data = data_fs_1, 
                      dist = "weibull")


plot_df <- summary(m2.fs1, type = "survival")

plot_df$`trans=12`$trans <- "instar 1 to 2"
plot_df$`trans=23`$trans <- "instar 2 to 3"
plot_df$`trans=34`$trans <- "instar 3 to 4"
plot_df$`trans=10`$trans <- "instar 1 to death"
plot_df$`trans=20`$trans <- "instar 2 to death"
plot_df$`trans=30`$trans <- "instar 3 to death"


plot_df2 <- rbind(plot_df$`trans=12`, plot_df$`trans=23`, plot_df$`trans=34`, plot_df$`trans=10`, plot_df$`trans=20`, plot_df$`trans=30`)

ggplot(data = plot_df2, aes(x = time, y = est)) +
  geom_line(aes(color = as.factor(trans))) +
  theme_bw() +
  labs(x = "Days since start of instar", y = "probability of remaining in current state")

tmat <- rbind(
  c(NA, 12, NA, NA, 10),
  c(NA, NA, 23, NA, 20),
  c(NA, NA, NA, 34, 30),
  c(NA, NA, NA, NA, 40),
  c(NA, NA, NA, NA, NA)
)

tgrid <- seq(0, 10, by  =1)

m2.pred <- msfit.flexsurvreg(m2.fs1, t = tgrid, trans = tmat)

m2.cox <- coxph(Surv(duration, status) ~ strata(trans), data = data_fs_1)
m3.cox <- coxph(Surv(start, end, status) ~ strata(trans), data = data_fs_1)




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
