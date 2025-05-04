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
library(lme4)
library(ggeffects)
library(sjPlot)

#### reading in data ####
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
            ordered = TRUE) %>% 
  filter(!(stage == "pupa" & rearing_year == 2024))

#### survivorship model -- logistic regression ####

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

plot_data <- surv_data %>% 
  group_by(stage, survival, after_june) %>% 
  summarize(count = n()) %>% 
  ungroup() %>% 
  group_by(stage, after_june) %>% 
  mutate(percent_surv = count / sum(count)) %>% 
  filter(survival == 1)

ggplot(plot_data, aes(x = stage, y = percent_surv)) +
         geom_col(aes(fill = after_june), alpha = 0.9, position = "dodge") +
  scale_size_continuous(limits = c(0.7, 1), range = c(0.02,20), breaks = c(0.75, 1.0)) +
  theme_bw()

surv_model <- glmer(survival ~ 
                      after_2018 +
                      after_june +
                      mean_temp_scaled * stage +
                      (1 | rearing_year),
                    data = surv_data,
                    family = binomial,
                    control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))
                    )
summary(surv_model)
ranef(surv_model)
resid(surv_model)

exp(fixef(surv_model))

AIC(surv_model)

plot_model(surv_model, type = "est", show.values = TRUE, value.offset = .3, vline.color = "grey70", vline.linetype = "dashed") + 
  theme_bw() +
  labs(title = NULL)

plot_model(surv_model_pupae, type = "est", show.values = TRUE, value.offset = .3, vline.color = "grey70", vline.linetype = "dashed") + 
  theme_bw() +
  labs(title = NULL)

pred_grid <- expand.grid(
  mean_temp_scaled = seq(-2.5, 2.5, by = 0.1),
  after_2018 = unique(surv_data$after_2018),
  after_june = unique(surv_data$after_june),
  stage = unique(surv_data$stage)
) %>% 
  mutate(mean_temp_unscaled = (mean_temp_scaled * mean_temp_sd) + mean_temp_mean)

pred_grid$predicted <- predict(surv_model, newdata = pred_grid, type = "response", re.form = NA)

# Plot predicted probabilities
ggplot(pred_grid %>% filter(after_2018 == "TRUE"), aes(x = mean_temp_unscaled, y = predicted, color = stage)) +
  geom_line(aes(linetype = after_june)) +
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = after_2018), alpha = 0.2, color = NA) +
  labs(x = "Mean Temperature (°C)",
       y = "Predicted Probability of Survival") +
  theme_bw()

# Generate predicted probabilities
preds <- ggpredict(surv_model, terms = c("after_2018", "after_june"))

emm <- emmeans(surv_model, ~ after_2018 + after_june, type = "response")
emm_df <- as.data.frame(emm)

ggplot(emm_df, aes(x = after_2018, y = prob)) +
  geom_col(position = "dodge") +
  labs(x = "After 2018", y = "Estimated Survival Probability") +
  theme_minimal()

ggplot(emm_df, aes(x = after_june, y = prob)) +
  geom_col(position = "dodge") +
  labs(x = "After June", y = "Estimated Survival Probability") +
  theme_minimal()

######## separate pupae model to make sure temperature relationship still holds (survivorship model) ########
surv_data_pupae <- temp_data %>% 
  mutate(survival = as.factor(if_else(stage2 == "death", 0, 1)),
         stage = factor(stage, ordered = F)) %>% 
  filter(stage == "pupa") %>% 
  mutate(
    rearing_year_scaled = scale(rearing_year),
    mean_temp_scaled = scale(mean_temp),
    after_2018 = as.factor(after_2018),
    after_june = as.factor(after_june))

# for unscaling later
mean_temp_sd_pupae <- sd(surv_data_pupae$mean_temp)
mean_temp_mean_pupae <- mean(surv_data_pupae$mean_temp)

surv_model_pupae <- glmer(survival ~ 
                            after_2018 +
                            after_june +
                            mean_temp_scaled +
                            (1 | rearing_year),
                          data = surv_data_pupae,
                          family = binomial,
                          control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))
)
summary(surv_model_pupae)



#### duration of stages + temperature model -- AFT ####
# TODO: try to make this multi state


multi_state_data <- temp_data %>%
  filter(!is.na(stage2)) %>%
  mutate(
    from_stage = stage,
    to_stage = stage2,
    status = if_else(stage2 == "death", 0, 1) # 1 = successfully transitioned, 0 = censored (death)
  )

keep_stages <- multi_state_data %>%
  group_by(from_stage) %>%
  summarize(deaths = sum(status == 0)) %>%
  filter(deaths > 5) %>%
  pull(from_stage)

multi_state_data_filtered <- multi_state_data %>%
  filter(from_stage %in% keep_stages)

multi_state_data_filtered %>%
  summarize(
    na_duration = sum(duration == 0),
    na_status = sum(status == 0),
    na_mean_temp = sum(mean_temp == 0)
  )

multi_state_model <- flexsurvreg(
  Surv(duration, status) ~ mean_temp + strata(from_stage),
  data = multi_state_data_filtered %>% filter(stage != "pupa"),
  dist = "weibull"
)
multi_state_model

# separate models for each stage
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
  mutate(stage = factor(stage, levels = c("First", "Second", "Third", "Fourth", "Fifth", "Pupa"), ordered = TRUE))

ggplot(res, aes(x = estimate, y = stage)) +
  geom_vline(xintercept = 1.00, col = "red", linetype = 2) +
  geom_segment(aes(x = conf.low, xend = conf.high), col = "darkgrey") +
  geom_point(size = 3, col = "darkgreen") +
  labs(x = "Hazard Ratio", y = "Stage") +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_model_fn <- function(model_to_plot, stage = "", temps = c(12, 14, 16), days = 15) {
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
    geom_line(aes(y = lcl), linetype = "dotted") +
    geom_line(aes(y = ucl), linetype = "dotted") +
    # geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = factor(mean_temp)), alpha = 0.4, linewidth = 0) +
    geom_line() +
    labs(x = paste0("Days after entering stage"),
         y = "Percent remaining in stage",
         color = "Mean temperature",
         title = stage) +
    theme_bw()
}

p1 <- plot_model_fn(model_first_instar, "First instar")
p2 <- plot_model_fn(model_second_instar, "Second instar")
p3 <- plot_model_fn(model_third_instar, "Third instar")
p4 <- plot_model_fn(model_fourth_instar, "Fourth instar")
p5 <- plot_model_fn(model_fifth_instar, "Fifth instar")

legend <- get_legend(p1)

plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"), 
  p3 + theme(legend.position = "none"), 
  p4 + theme(legend.position = "none"),
  p5 + theme(legend.position = "none"),
  legend)

pp <- plot_model_fn(models_surv$pupa, "Pupa", temps = c(9, 10, 11, 12), days = 450)

# raw data plots
r1 <- ggplot(data_first_instar, aes(x = duration, y = mean_temp)) +
  geom_jitter(alpha = 0.3, aes(colour = as.factor(rearing_year))) +
  geom_smooth(method = "lm") +
  labs(title = "First Instar") +
  ylim(9, 25) +
  theme_bw()

r2 <- ggplot(data_second_instar %>% filter(status == 1), aes(x = duration, y = mean_temp)) +
  geom_jitter(alpha = 0.3, aes(colour = as.factor(rearing_year))) +
  geom_smooth(method = "lm") +
  labs(title = "Second Instar") +
  ylim(9, 25) +
  theme_bw()

r3 <- ggplot(data_third_instar %>% filter(status == 1), aes(x = duration, y = mean_temp)) +
  geom_jitter(alpha = 0.3, aes(colour = as.factor(rearing_year))) +
  geom_smooth(method = "lm") +
  labs(title = "Third Instar") +
  ylim(9, 25) +
  theme_bw()

r4 <- ggplot(data_fourth_instar %>% filter(status == 1), aes(x = duration, y = mean_temp)) +
  geom_jitter(alpha = 0.3, aes(colour = as.factor(rearing_year))) +
  geom_smooth(method = "lm") +
  labs(title = "Fourth Instar") +
  ylim(9, 25) +
  theme_bw()

r5 <- ggplot(data_fifth_instar %>% filter(status == 1), aes(x = duration, y = mean_temp)) +
  geom_jitter(alpha = 0.3, aes(colour = as.factor(rearing_year))) +
  geom_smooth(method = "lm") +
  labs(title = "Fifth Instar") +
  ylim(9, 25) +
  theme_bw()

rp <- ggplot(data_pupa %>% filter(status == 1), aes(x = duration, y = mean_temp)) +
  geom_jitter(alpha = 0.3, aes(colour = as.factor(rearing_year))) +
  geom_smooth(method = "lm") +
  labs(title = "Pupa") +
  ylim(9, 15) +
  theme_bw()

legend <- get_legend(
  r1 + labs(color = "Rearing Year")
  )

rps <- plot_grid(r1 + theme(legend.position = "none"), 
          r2 + theme(legend.position = "none"), 
          r3 + theme(legend.position = "none"), 
          r4 + theme(legend.position = "none"), 
          r5 + theme(legend.position = "none"), 
          rp + theme(legend.position = "none"))

plot_grid(rps, legend, nrow = 1, rel_widths = c(3, 1))

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
  xlim(9, 15) +
  theme_bw()


plot_grid(t1 + theme(legend.position = "none"), 
          t2 + theme(legend.position = "none") + labs(y = NULL), 
          t3 + theme(legend.position = "none") + labs(y = NULL), 
          t4 + theme(legend.position = "none"), 
          t5 + theme(legend.position = "none") + labs(y = NULL), 
          tp + theme(legend.position = "none") + labs(y = NULL))


#### archive -- survival models ####

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

# fit_year <- survfit(Surv(start, end, stage2) ~ rearing_year, data = larvae_data, id = imb_id)


fit <- fit_null # put desired fit here to plot
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