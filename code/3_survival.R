library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(ggplot2)
library(survival)
library(flexsurv)
library(cowplot)
library(lme4)
library(sjPlot)
library(stringr)
library(RColorBrewer)

#### reading in data ####
temp_data <- read_csv(here::here("data", "mean_temp_data.csv")) %>% 
  mutate(
    imb_id = as.factor(imb_id),
    host_plant = as.factor(host_plant),
    rearing_year = as.integer(rearing_year),
    stage = as.factor(stage)
  ) %>% 
  mutate_at(vars(contains("stage")), 
            factor,
            levels = c("first", "second", "third", "fourth", "fifth", "pupa", "adult", "death"),
            ordered = TRUE) %>% 
  filter(!(stage == "pupa" & rearing_year == 2024),
         !(stage == "pupa" & duration ==1))

# summary data table
summary <- temp_data %>% filter(stage2 != "death") %>% 
  group_by(stage) %>% 
  summarise(mean_duration = mean(duration),
            sd_duration = sd(duration))
  
# formatting data for model
surv_data <- temp_data %>% 
mutate(survival = as.factor(if_else(stage2 == "death", 0, 1)),
       stage = factor(stage, ordered = F),
       rearing_year_scaled = scale(rearing_year),
       mean_temp_scaled = scale(mean_temp),
       after_2018 = as.factor(after_2018),
       after_june = as.factor(after_june))

# for unscaling later
mean_temp_sd <- sd(surv_data$mean_temp)
mean_temp_mean <- mean(surv_data$mean_temp)

#### survfit plots for visualization of stages through time ####

fit_null <- survfit(Surv(start, end, stage2) ~ 1, data = surv_data, id = imb_id)
fit_larvae <- survfit(Surv(start, end, stage2) ~ 1, data = surv_data %>% filter(stage != "pupa"), id = imb_id)

fit <- fit_larvae # put desired fit here to plot
# fit <- fit_null

# plot(fit, col = palette, lwd = 2, noplot = NULL,
#      xlab = "Days after entering the first instar",
#      ylab = "Probability in state")
# legend(35, .75, c("1st instar", "2nd instar", "3rd instar", "4th instar", "5th instar", "pupa", "death"), col = c(1, 2, 3, 4, 5, 6, 7), lty = 1)

df <- data.frame(
  time = fit$time,
  prob = fit$pstate
) %>% 
  select(-prob.adult) %>% 
  rename(
    prob.first = prob..s0.
  ) %>% 
  pivot_longer(
    cols = starts_with("prob"),
    names_to = "stage",
    values_to = "probability"
  ) %>% 
  mutate(
    stage = fct_inorder(str_sub(stage, start = 6))
  )

ggplot(df, aes(x = time, y = probability, color = stage)) +
  geom_smooth(se = F, span = 0.1) +
  lims(y = c(0, 1.1)) +
  labs(x = "Days after entering the first instar",
       y = "Probability in state") +
  theme_bw() +
  scale_color_manual(
    values = brewer.pal(8, "Dark2"), #set1
    name = "Developmental Stage",
    labels = c("First Instar", "Second Instar", "Third Instar", "Fourth Instar", "Fifth Instar", "Pupa", "Death"))

#### survivorship model -- logistic regression ####

# glmer model
surv_model <- glmer(survival ~ 
                      after_2018 +
                      after_june +
                      stage +
                      (1 | rearing_year),
                    data = surv_data,
                    family = binomial,
                    control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))
                    )
summary(surv_model)

# forest plot of odds ratios by stage, after 2018, and after june (using sjPlot)
plot_model(surv_model, type = "est", show.values = TRUE, value.offset = .3, vline.color = "grey70", vline.linetype = "dashed") + 
  theme_bw() +
  labs(title = NULL)

# plots showing likelihood of survival through each stage depending on collection month and year
pred_grid <- expand.grid(
  after_2018 = unique(surv_data$after_2018),
  after_june = unique(surv_data$after_june),
  stage = unique(surv_data$stage)
)

predictions <- predict(surv_model, newdata = pred_grid, re.form = NA, type = "response")


# this takes forever to run, I saved it to boot_glmer_data.csv

# predict_fun <- function(model) {
#   predict(model, newdata = pred_grid, re.form = NA, type = "response")
# }
# 
# boot_results <- bootMer(
#   surv_model, 
#   predict_fun, 
#   nsim = 500, 
#   use.u = FALSE, 
#   type = "parametric",
#   verbose = TRUE,
#   parallel = "multicore",
#   ncpus = 4
# )

# res <- data.frame(pred_grid, boot_results$t0, confint(boot_results)) %>% 
#   rename(est = boot_results.t0,
#          lower = X2.5..,
#          upper = X97.5..)

# write_csv(res, here::here("data", "boot_glmer_data.csv"))


# run this instead to read the data
res <- read_csv(here::here("data", "boot_glmer_data.csv")) %>% 
  mutate(
    stage = fct_inorder(stage)
  )

plot_data <- res %>% 
  filter(after_2018 == "TRUE", after_june == "FALSE")

# plot_data <- res %>% 
#   filter(after_2018 == "FALSE", after_june == "TRUE")

stage_plot <- ggplot(plot_data, aes(x = stage, y = est)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  # geom_crossbar(aes(ymin = lower, ymax = upper), fill = "grey") +
  labs(x = "Developmental Stage",
       y = "Predicted Probability of Survival") +
  lims(y = c(0.8, 1)) +
  theme_bw() +
  scale_x_discrete(labels = c("First", "Second", "Third", "Fourth", "Fifth", "Pupa")) +
  theme(text = element_text(size = 14))
stage_plot

cumulative_surv_data <- res %>% 
  group_by(after_2018, after_june) %>% 
  mutate(
    cum_surv = cumprod(est),
    cum_lower = cumprod(lower),
    cum_upper = cumprod(upper),
    stage_num = as.numeric(stage)
  ) %>% 
  ungroup() %>% 
  filter(after_2018 == TRUE, after_june == FALSE)
  
cum_surv_plot <- ggplot(cumulative_surv_data, aes(x = stage_num, y = cum_surv)) +
  geom_ribbon(aes(ymin = cum_lower, ymax = cum_upper), alpha = 0.6, fill = "#817E9F") +
  geom_line(linewidth = 1.5) +
  labs(x = "Developmental Stage",
       y = "Predicted Cumulative Probability of Survival") +
  lims(y = c(0.80, 1)) +
  scale_x_continuous(breaks = 1:6, labels = c("First", "Second", "Third", "Fourth", "Fifth", "Pupa")) +
  theme_bw() +
  theme(text = element_text(size = 14))

cum_surv_plot

cowplot::plot_grid(stage_plot, cum_surv_plot)

# temperature and duration EDA plots -- these go with AFT model

surv_data_pupae <- surv_data %>% 
  filter(stage == "pupa")

pd2 <- surv_data_pupae %>% filter(stage2 == "adult", rearing_year != 2023) %>% mutate(t2 = as.ordered(round(mean_temp, 1)))

ggplot(surv_data_pupae, aes(x = duration)) +
  geom_bar(aes(fill = stage2), alpha = 0.6) +
  theme_bw()


ggplot(pd2, aes(x = duration, y = mean_temp)) +
  geom_point(aes(color = factor(rearing_year)), alpha = 0.5) +
  geom_smooth(aes(color = factor(rearing_year)), se = F) +
  theme_bw() +
  labs(x = "Duration of pupal stage",
       y = "Mean temperature",
       color = "Rearing year")

ggplot(temp_data %>% filter(stage != "pupa", stage2 != "death"), aes(x = duration, y = mean_temp)) +
  geom_jitter(aes(color = stage), alpha = 0.1, size = 2) +
  geom_smooth(aes(color = stage), method = "lm", se = F) +
  theme_bw() +
  labs(x = "Duration of instar stage",
       y = "Mean temperature",
       color = "Larval instar")


#### duration of stages + temperature model -- AFT ####
subset_data <- function(stage_name) {
  temp_data %>% 
    group_by(imb_id) %>% 
    filter(stage == stage_name) %>% 
    mutate(status = if_else(stage2 == "death", 0, 1))
}

# duration models
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
  
res <- rbind(m1, m2, m3, m4, m5, mp) %>% 
  filter(term == "mean_temp") %>% 
  mutate(stage = factor(stage, levels = c("First", "Second", "Third", "Fourth", "Fifth", "Pupa"), ordered = TRUE),
         stage = fct_rev(stage))

# AFT results plot showing how much 1 degree of temp increase multiplies the length of time in the stage
ggplot(res, aes(x = estimate, y = stage)) +
  geom_vline(xintercept = 1.00, col = "#BD210F", linetype = 2, linewidth = 0.6) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), col = "#817E9F") +
  geom_point(size = 3, col = "black") +
  labs(x = "Stage", y = "Time Ratio") +
  theme_bw(base_size = 15)

palette <- c("#319CD2", "#474959", "#C2190A", "black")

# plots showing likelihood of remaining in stage over time for different values of temperature
plot_model_fn <- function(model_to_plot, stage = "", temps = c(11, 14, 17), days = 15, ylab = NULL) {
  # generating df for plotting
  newdata <- data.frame(mean_temp = temps)
  time_range <- seq(0, days, by = 1)
  
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
    # geom_line(aes(y = lcl), linetype = "dotted") +
    # geom_line(aes(y = ucl), linetype = "dotted") +
    geom_line(size = 1.3) +
    scale_color_manual(
      values = palette,
      name = "Mean Temperature",
      labels = c("11 °C", "14 °C", "17 °C")) +
    labs(x = "Days after entering stage",
         y = ylab,
         title = stage) +
    theme_bw() +
    theme(text = element_text(size = 14))
}

p1 <- plot_model_fn(model_first_instar, "First instar", ylab = "Probability of remaining in stage")
p2 <- plot_model_fn(model_second_instar, "Second instar")
p3 <- plot_model_fn(model_third_instar, "Third instar")
p4 <- plot_model_fn(model_fourth_instar, "Fourth instar", ylab = "Probability of remaining in stage")
p5 <- plot_model_fn(model_fifth_instar, "Fifth instar")

legend <- get_legend(p1)

cowplot::plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"), 
  p3 + theme(legend.position = "none"), 
  p4 + theme(legend.position = "none"),
  p5 + theme(legend.position = "none"),
  legend)

pp <- plot_model_fn(model_pupa, "Pupa", temps = c(9, 10, 11, 12), days = 450)
pp
