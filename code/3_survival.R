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
  filter(!(stage == "pupa" & rearing_year == 2024),
         !(stage == "pupa" & duration ==1))

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

# plot_data <- surv_data %>% 
#   group_by(stage, survival, after_june) %>% 
#   summarize(count = n()) %>% 
#   ungroup() %>% 
#   group_by(stage, after_june) %>% 
#   mutate(percent_surv = count / sum(count)) %>% 
#   filter(survival == 1)

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
ranef(surv_model)
resid(surv_model)

exp(fixef(surv_model))

AIC(surv_model)

plot_model(surv_model, type = "est", show.values = TRUE, value.offset = .3, vline.color = "grey70", vline.linetype = "dashed") + 
  theme_bw() +
  labs(title = NULL)

pred_grid <- expand.grid(
  after_2018 = unique(surv_data$after_2018),
  after_june = unique(surv_data$after_june),
  stage = unique(surv_data$stage)
)
pred_grid$predicted <- predict(surv_model, newdata = pred_grid, type = "response", re.form = NA)

plot_data <- pred_grid %>% 
  filter(after_2018 == "TRUE") %>% 
  arrange(after_june) %>% 
  mutate(after_june = if_else(after_june == "TRUE", "During/After June", "Before June"))

# Plot predicted probabilities
ggplot(plot_data, aes(x = stage, y = predicted, group = after_june)) +
  geom_point(aes(color = after_june)) +
  geom_segment(aes(xend = if_else(stage == "pupa", stage, lead(stage)), yend = if_else(stage == "pupa", predicted, lead(predicted)), color = after_june)) +
  labs(x = "Developmental Stage",
       y = "Predicted Probability of Survival",
       color = "Collection Date") +
  theme_bw()

plot_data <- pred_grid %>% 
  filter(after_june == "FALSE") %>% 
  arrange(after_2018) %>% 
  mutate(after_2018 = if_else(after_2018 == "TRUE", "During/After 2018", "Before 2018"))

# Plot predicted probabilities
ggplot(plot_data, aes(x = stage, y = predicted, group = after_2018)) +
  geom_point(aes(color = after_2018)) +
  geom_segment(aes(xend = if_else(stage == "pupa", stage, lead(stage)), yend = if_else(stage == "pupa", predicted, lead(predicted)), color = after_2018)) +
  labs(x = "Developmental Stage",
       y = "Predicted Probability of Survival",
       color = "Collection Date") +
  theme_bw()

ggplot(surv_data_pupae, aes(x = duration)) +
  geom_bar(aes(fill = stage2), alpha = 0.6) +
  theme_bw()

pd2 <- surv_data_pupae %>% filter(stage2 == "adult", rearing_year != 2023) %>% mutate(t2 = as.ordered(round(mean_temp, 1)))

ggplot(pd2, aes(x = duration)) +
  geom_bar(aes(fill = t2)) +
  facet_wrap(~rearing_year) +
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
# TODO: try to make this multi state
library(mstate)

states <- c("first", "second", "third", "fourth", "fifth", "pupa", "death")

tmat <- transMat(
  list(
    c(2, 7),    # from first
    c(3, 7),     # from second
    c(4, 7),    # from third
    c(5),     # from fourth --> only to 5th instar to make model converge (hopefully lol)
    c(6, 7),      # from fifth
    c(),              # from pupa
    c()                      # death is absorbing
  ),
  names = states
)

tmat

# Rename to match `msprep()` naming convention
msdata <- multi_state_data %>%
  filter(stage %in% states, stage2 %in% states) %>%
  mutate(
    from = as.character(stage),
    to = as.character(stage2),
    id = imb_id
  )

# Recode from and to as factors with correct levels
msdata$from <- factor(msdata$from, levels = states)
msdata$to   <- factor(msdata$to,   levels = states)

msdata$trans <- tmat[cbind(as.numeric(msdata$from), as.numeric(msdata$to))]

msdata <- msdata %>% 
  filter(!is.na(trans))

fit <- flexsurvreg(
  Surv(start, end, status) ~ mean_temp + factor(stage, ordered = F),
  data = msdata,
  dist = "weibull"
)
fit
plot(fit)

# Define temperatures at which you want to plot survival curves
temps <- c(10, 15, 20)  # example temperatures in degrees

# Define stages to plot (all levels except baseline)
stages <- levels(factor(msdata$stage))

# Create a new data frame for predictions with all combinations of temp and stage
newdata <- expand.grid(
  mean_temp = temps,
  stage = stages
)

# Predict the survival curves for each row in newdata
# Use the flexsurv package's summary.flexsurvreg with type = "survival"
# and specify a sequence of times for the survival curve

time_seq <- seq(0, 50, by = 0.5)  # adjust max time based on your data

# For each row in newdata, get survival probabilities over time
surv_list <- lapply(1:nrow(newdata), function(i) {
  summary(fit, newdata = newdata[i, ], type = "survival", t = time_seq)[[1]] %>%
    as.data.frame() %>%
    mutate(
      time = time_seq,
      mean_temp = newdata$mean_temp[i],
      stage = newdata$stage[i]
    )
})

# Combine all survival data frames
surv_df <- bind_rows(surv_list)

# Plot using ggplot2
ggplot(surv_df, aes(x = time, y = est, color = factor(stage), linetype = factor(mean_temp))) +
  geom_line(size = 1) +
  labs(
    x = "Time",
    y = "Survival Probability",
    color = "Stage",
    linetype = "Mean Temperature (°C)",
    title = "Predicted Survival Curves by Stage and Temperature"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right"
  )


#### duration models that actually work -- one for each stage ####
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
    geom_line(aes(y = lcl), linetype = "dotted") +
    geom_line(aes(y = ucl), linetype = "dotted") +
    # geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = factor(mean_temp)), alpha = 0.4, linewidth = 0) +
    geom_line() +
    labs(x = paste0("Days after entering stage"),
         y = ylab,
         color = "Mean temperature",
         title = stage) +
    theme_bw()
}

p1 <- plot_model_fn(model_first_instar, "First instar", ylab = "Percent remaining in stage")
p2 <- plot_model_fn(model_second_instar, "Second instar")
p3 <- plot_model_fn(model_third_instar, "Third instar")
p4 <- plot_model_fn(model_fourth_instar, "Fourth instar", ylab = "Percent remaining in stage")
p5 <- plot_model_fn(model_fifth_instar, "Fifth instar")

legend <- get_legend(p1)

cowplot::plot_grid(
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


fit <- fit_larvae # put desired fit here to plot
print(fit)

fit$transitions

plot(fit, col = c(1, 2, 3, 4, 5, 6, 1, 7), noplot = NULL,
     xlab = "Days after entering the first instar",
     ylab = "Probability in state")
legend(35, .75, c("1st instar", "2nd instar", "3rd instar", "4th instar", "5th instar", "pupa", "death"), col = c(1, 2, 3, 4, 5, 6, 7), lty = 1)


cox_dd <- coxph(Surv(start_dd, cum_degree_days, stage2) ~ rearing_year, data = data, id = imb_id)

cox_time <- coxph(Surv(start, end, stage2) ~ rearing_year, data = data, id = imb_id, model = T, x = T, y = T)
cox_time$transitions

summary(cox_time)
summary(cox_dd)



cox.zph(cox_time, transform = "identity")
cox.zph(cox_dd, transform = "identity")