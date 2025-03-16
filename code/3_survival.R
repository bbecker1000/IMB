library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)
library(ggplot2)
library(survival)
library(icenReg)
library(flexsurv)
library(cowplot)
library(purrr)


data <- read_csv(here::here("data", "dd_data.csv")) %>% 
  mutate(
    imb_id = as.factor(imb_id),
    host_plant = as.factor(host_plant),
    # overall_survival = as.factor(overall_survival),
    rearing_year = as.integer(rearing_year),
    stage = as.factor(stage)
  ) %>% 
  mutate_at(vars(contains("stage")), 
            factor,
            levels = c("first", "second", "third", "fourth", "fifth", "pupa", "adult", "death"),
            ordered = TRUE) 

larvae_data <- data %>% 
  group_by(imb_id) %>% 
  mutate(
    survives_to_pupa = "pupa" %in% stage2,
    to_exclude = if_else(survives_to_pupa & stage2 == 'death', TRUE, FALSE)
  ) %>% 
  filter(!to_exclude,
         # stage2 != "pupa", 
         stage2 != "adult")


# new approach -- one model for each transition state with mean temperature and year as covariates

temp_data <- read_csv(here::here("data", "mean_temp_data.csv")) %>% 
  mutate(
    imb_id = as.factor(imb_id),
    host_plant = as.factor(host_plant),
    # overall_survival = as.factor(overall_survival),
    rearing_year = as.integer(rearing_year),
    stage = as.factor(stage)
  ) %>% 
  mutate_at(vars(contains("stage")), 
            factor,
            levels = c("first", "second", "third", "fourth", "fifth", "pupa", "adult", "death"),
            ordered = TRUE) 

data_first_instar <- temp_data %>% 
  group_by(imb_id) %>% 
  filter(stage == "first" | (
    stage == "death" & lag(stage) == "first"
  )) %>% 
  mutate(status = if_else(stage == "first", 1, 0))

fit_first_instar <- survfit(Surv(start, end, status) ~ 1, data = data_first_instar)
print(fit_first_instar)
plot(fit_first_instar,
     xlab = "days after entering first instar",
     ylab = "percent in first instar")

model_first_instar <- coxph(Surv(start, end, status) ~ mean_temp, data = data_first_instar)
summary(model_first_instar)

z <- cox.zph(model_first_instar, transform = "identity")
print(z)
plot(z)

model_first_instar <- flexsurvreg(Surv(duration, status) ~ mean_temp, data = data_first_instar, dist = "weibull")
summary(model_first_instar)

# generating df for plotting
newdata <- data.frame(mean_temp = c(10, 15, 20))
time_range <- seq(0, 15, by = 1)

pred_list <- summary(
  model_first_instar,
  newdata = newdata,
  type = "survival",
  t = time_range,
  ci = TRUE
)

pred_df <- map2_dfr(
  pred_list, 
  newdata$mean_temp, 
  ~ data.frame(
    time       = .x$time,
    surv       = .x$est,
    lcl = .x$lcl,
    ucl = .x$ucl,
    mean_temp  = .y
  )
)

ggplot(pred_df, aes(x = time, y = surv, color = factor(mean_temp))) +
  geom_line(aes(y = lcl), linetype = "dotted") +
  geom_line(aes(y = ucl), linetype = "dotted") +
  # geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = factor(mean_temp)), alpha = 0.4, linewidth = 0) +
  geom_line() +
  labs(x = "Days after entering first instar",
       y = "Percent remaining in first instar",
       color = "Mean temperature") +
  theme_bw()

plot(model_first_instar, 
     newdata = newdata,
     col = c(1, 2, 3),
     xlab = "Days after entering first instar",
     ylab = "Percent in first instar")

# gonna try to just use the survival package because it also supports multistate models
# if using the survival package, groups (i.e. RHS of the model) may not be time dependent
# so we must use degree days as the time variable

# TODO: to make this work, i need starting degree days and ending degree days
# for now I will just use a null model with time, and if it works I'll modify degree days

fit_pre2018 <- survfit(Surv(start, end, stage2) ~ 1, data = data %>% filter(!after_2018), id = imb_id)
fit_post2018 <- survfit(Surv(start, end, stage2) ~ 1, data = data %>% filter(after_2018), id = imb_id)

fit_null <- survfit(Surv(start, end, stage2) ~ 1, data = data, id = imb_id)

fit_null_dd <- survfit(Surv(start_dd, cum_degree_days, stage2) ~ 1, data = data, id = imb_id)

fit_2018 <- survfit(Surv(start_dd, cum_degree_days, stage2) ~ after_2018, data = data, id = imb_id)

fit_larvae <- survfit(Surv(start, end, stage2) ~ 1, data = larvae_data, id = imb_id)

fit_larvae_dd <- survfit(Surv(start_dd, cum_degree_days, stage2) ~ 1, data = larvae_data, id = imb_id)

fit_larvae_2018 <- survfit(Surv(start_dd, cum_degree_days, stage2) ~ after_2018, data = larvae_data, id = imb_id)

fit_year <- survfit(Surv(start, end, stage2) ~ rearing_year, data = larvae_data, id = imb_id)


fit <- fit_year # put desired fit here to plot
print(fit)

fit$transitions

plot(fit, col = c(1, 2, 3, 4, 5, 6, 7), noplot = NULL,
     xlab = "degree days after entering the first instar",
     ylab = "probability in state")
legend(150, .75, c("1st instar", "2nd instar", "3rd instar", "4th instar", "5th instar", "pupa", "death"), col = c(1, 3, 4, 5, 6, 7, 2), lty = 1)


cox_dd <- coxph(Surv(start_dd, cum_degree_days, stage2) ~ rearing_year, data = data, id = imb_id)

cox_time <- coxph(Surv(start, end, stage2) ~ rearing_year, data = data, id = imb_id, model = T, x = T, y = T)
cox_time$transitions

summary(cox_time)
summary(cox_dd)



cox.zph(cox_time, transform = "identity")
cox.zph(cox_dd, transform = "identity")

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


