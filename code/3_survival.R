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
  filter(!to_exclude, stage2 != "adult")


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

subset_data <- function(stage_name) {
  temp_data %>% 
    group_by(imb_id) %>% 
    filter(stage == stage_name) %>% 
    mutate(status = if_else(stage2 == "death", 0, 1))
}

data_first_instar <- subset_data(stage_name = "first")
data_second_instar <- subset_data("second")
data_third_instar <- subset_data("third")
data_fourth_instar <- subset_data("fourth")
data_fifth_instar <- subset_data("fifth")
data_pupa <- subset_data("pupa")

model_first_instar <- flexsurvreg(Surv(duration, status) ~ mean_temp, data = data_first_instar, dist = "weibull")
m1 <- tidy(model_first_instar, conf.int = TRUE, transform = "coefs.exp") %>% mutate(stage = "First")

model_second_instar <- flexsurvreg(Surv(duration, status) ~ mean_temp, data = data_second_instar, dist = "weibull")
m2 <- tidy(model_second_instar, conf.int = TRUE, transform = "coefs.exp") %>% mutate(stage = "Second")

model_third_instar <- flexsurvreg(Surv(duration, status) ~ mean_temp, data = data_third_instar, dist = "weibull")
m3 <- tidy(model_third_instar, conf.int = TRUE, transform = "coefs.exp") %>% mutate(stage = "Third")

model_fourth_instar <- flexsurvreg(Surv(duration, status) ~ mean_temp, data = data_fourth_instar, dist = "weibull")
m4 <- tidy(model_fourth_instar, conf.int = TRUE, transform = "coefs.exp") %>% mutate(stage = "Fourth")

model_fifth_instar <- flexsurvreg(Surv(duration, status) ~ mean_temp, data = data_fifth_instar, dist = "weibull")
m5 <- tidy(model_fifth_instar, conf.int = TRUE, transform = "coefs.exp") %>% mutate(stage = "Fifth")

model_pupa <- flexsurvreg(Surv(duration, status) ~ mean_temp, data = data_pupa, dist = "weibull")
mp <- tidy(model_pupa, conf.int = TRUE, transform = "coefs.exp") %>% mutate(stage = "Pupa")

model_res <- rbind(m1, m2, m3, m4, m5, mp) %>% 
  filter(term == "mean_temp") %>% 
  mutate(stage = factor(stage, levels = c("First", "Second", "Third", "Fourth", "Fifth", "Pupa"), ordered = TRUE))

ggplot(model_res, aes(x = estimate, y = stage)) +
  geom_vline(xintercept = 1.00, col = "red", linetype = 2) +
  geom_segment(aes(x = conf.low, xend = conf.high), col = "darkgrey") +
  geom_point(size = 3, col = "darkgreen") +
  labs(x = "Hazard Ratio", y = "Stage") +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_model_fn <- function(model_to_plot, stage = "", temps = c(12, 14, 16)) {
  # generating df for plotting
  newdata <- data.frame(mean_temp = temps)
  time_range <- seq(0, 15, by = 1)
  
  pred_list <- summary(
    model_to_plot,
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
    labs(x = paste0("Days after entering ", stage, " instar"),
         y = "Percent remaining in instar",
         color = "Mean temperature") +
    theme_bw()
  }

p1 <- plot_model_fn(model_first_instar, "first")
p2 <- plot_model_fn(model_second_instar, "second")
p3 <- plot_model_fn(model_third_instar, "third")
p4 <- plot_model_fn(model_fourth_instar, "fourth")
p5 <- plot_model_fn(model_fifth_instar, "fifth")

legend <- get_legend(p1)

plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"), 
  p3 + theme(legend.position = "none"), 
  p4 + theme(legend.position = "none"), 
  p5 + theme(legend.position = "none"),
  legend)

# comparisons of mean temps for surviving vs non surviving individuals in each stage
t1 <- ggplot(data_first_instar, aes(x = mean_temp, fill = stage2)) +
  geom_density(alpha = 0.4) +
  labs(title = "First instar") +
  xlim(9, 25) +
  theme_bw()

t2 <- ggplot(data_second_instar, aes(x = mean_temp, fill = stage2)) +
  geom_density(alpha = 0.4) +  
  labs(title = "Second instar") +
  xlim(9, 25) +
  theme_bw()

t3 <- ggplot(data_third_instar, aes(x = mean_temp, fill = stage2)) +
  geom_density(alpha = 0.4) +
  labs(title = "Third instar") +
  xlim(9, 25) +
  theme_bw()

t4 <- ggplot(data_fourth_instar, aes(x = mean_temp, fill = stage2)) +
  geom_density(alpha = 0.4) +
  labs(title = "Fourth instar") +
  xlim(9, 25) +
  theme_bw()

t5 <- ggplot(data_fifth_instar, aes(x = mean_temp, fill = stage2)) +
  geom_density(alpha = 0.4) +
  labs(title = "Fifth instar") +
  xlim(9, 25) +
  theme_bw()

tp <- ggplot(data_pupa, aes(x = mean_temp, fill = stage2)) +
  geom_density(alpha = 0.4) +
  labs(title = "Pupa") +
  xlim(9, 25) +
  theme_bw()


plot_grid(t1 + theme(legend.position = "none"), 
          t2 + theme(legend.position = "none") + labs(y = NULL), 
          t3 + theme(legend.position = "none") + labs(y = NULL), 
          t4 + theme(legend.position = "none"), 
          t5 + theme(legend.position = "none") + labs(y = NULL), 
          tp + theme(legend.position = "none") + labs(y = NULL))

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


