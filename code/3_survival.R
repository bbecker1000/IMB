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
library(marginaleffects)

# file paths
DATA_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)


##############################
# 1) Data import & cleaning
##############################
temp_data <- read_csv(file.path(DATA_DIR, "mean_temp_data.csv")) %>% 
  mutate(
    imb_id = as.factor(imb_id),
    host_plant = as.factor(host_plant),
    rearing_year = as.integer(rearing_year),
    stage = factor(stage)
  ) %>% 
  mutate_at(vars(contains("stage")), 
            factor,
            levels = c("first", "second", "third", "fourth", "fifth", "pupa", "adult", "death"),
            ordered = TRUE) %>% 
  # excluding non-eclosed pupae and mispupations
  filter(!(stage == "pupa" & rearing_year == 2024),
         !(stage == "pupa" & duration ==1))

# Quick summary table (exclude terminal death rows for duration stats)
stage_duration_summary <- temp_data %>%
  filter(stage2 != "death") %>%
  group_by(stage) %>%
  summarise(
    n = dplyr::n(),
    mean_duration = mean(duration, na.rm = TRUE),
    sd_duration = sd(duration, na.rm = TRUE),
    .groups = "drop"
  )




##############################
# 2) Survfit Graph: Creating data frame for generating figure
##############################

# formatting data for model
surv_data <- temp_data %>%
  mutate(
    # Binary survival indicator for GLMM: 1 = survived this stage (i.e., not death next), 0 = died
    survival = as.integer(if_else(stage2 == "death", 0, 1)),
    # Unordered factor for fixed effects where order does not imply linear trend
    stage_unordered = factor(as.character(stage), levels = c("first","second","third","fourth","fifth","pupa","adult","death"), ordered = FALSE),
    # Scaled covariates (as numeric vectors)
    rearing_year_scaled = as.numeric(scale(rearing_year)),
    mean_temp_scaled = as.numeric(scale(mean_temp)),
    # Ensure logicals are factors where needed by modeling conveniences
    after_2018 = factor(after_2018),
    after_june = factor(after_june)
  )

# Save means/SDs for later back-transform reference
mean_temp_mean <- mean(surv_data$mean_temp, na.rm = TRUE)
mean_temp_sd <- sd(surv_data$mean_temp, na.rm = TRUE)


##############################
# 3) Multi-state visualization with survfit
##############################
# Visualizes occupancy probability over time across states (pstate),

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


##############################
# 4) Survivorship GLMM (stage-level survival yes/no)
##############################
# Mixed-effects logistic regression: survival ~ after_2018 + after_june + stage + (1 | rearing_year)

surv_model <- glmer(survival ~ 
                      after_2018 +
                      after_june +
                      stage_unordered +
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
  stage_unordered = unique(surv_data$stage_unordered)
)

pred_grid$pred <- predict(surv_model, newdata = pred_grid, re.form = NA, type = "response")


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

# these plots but with marginaleffects


##############################
# 5) EDA — temperature vs duration
##############################

surv_data_pupae <- surv_data %>% 
  filter(stage == "pupa")

p_pupa_hist <- ggplot(surv_data_pupae, aes(x = duration, fill = stage2)) +
  geom_bar(aes(fill = stage2), alpha = 0.6) +
  theme_bw() + labs(x = "Duration of pupal stage", y = "Count", fill = "Next state")

p_pupa_scatter <- surv_data %>%
  filter(stage_unordered == "pupa", stage2 == "adult", rearing_year != 2023) %>%
  ggplot(aes(x = duration, y = mean_temp, color = factor(rearing_year))) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = FALSE) +
  theme_bw() + labs(x = "Duration of pupal stage", y = "Mean temperature", color = "Rearing year")

# Larval instars EDA
p_instars <- temp_data %>%
  filter(stage != "pupa", stage2 != "death") %>%
  ggplot(aes(x = duration, y = mean_temp, color = stage)) +
  geom_jitter(alpha = 0.1, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() + labs(x = "Duration of instar stage", y = "Mean temperature", color = "Larval instar")

##############################
# 6) AFT (Weibull) models: duration ~ mean_temp by stage
##############################
# We model each stage separately with flexsurvreg and extract time ratios (exp of coefficients)

subset_data <- function(stage_name) {
  temp_data %>% 
    group_by(imb_id) %>% 
    filter(stage == stage_name) %>% 
    mutate(status = if_else(stage2 == "death", 0, 1))
}

fit_aft_by_stage <- function(stage_name) {
  dat <- subset_data(stage_name)
  if (nrow(dat) < 5) return(NULL)
  model <- flexsurvreg(Surv(duration, status) ~ mean_temp, data = dat, dist = "weibull")
  broom::tidy(model, conf.int = TRUE, transform = "coefs.exp") %>%
    mutate(stage_label = str_to_title(stage_name))
}

stages_vec <- c("first","second","third","fourth","fifth","pupa")

aft_tidy <- map_dfr(stages_vec, fit_aft_by_stage) %>%
  filter(term == "mean_temp") %>%
  mutate(stage_label = factor(stage_label, levels = c("First","Second","Third","Fourth","Fifth","Pupa"), ordered = TRUE)) %>%
  arrange(stage_label)

# time ratio by stage plot
p_aft_tr <- ggplot(aft_tidy, aes(x = estimate, y = forcats::fct_rev(stage_label))) +
  geom_vline(xintercept = 1.00, color = "#BD210F", linetype = 2, linewidth = 0.6) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), color = "#817E9F") +
  geom_point(size = 3, color = "black") +
  labs(y = "Stage", x = "Time Ratio (per +1 °C)") +
  theme_bw(base_size = 15)
# ggsave(file.path(FIG_DIR, "aft_time_ratios.png"), p_aft_tr, width = 7, height = 5, dpi = 300)

# Survival curves from AFT models by temperature values
PLOT_PALETTE <- c("#319CD2", "#474959", "#C2190A", "black")

plot_aft_survival_curves <- function(model, stage_title = "", temps = c(11, 14, 17), days = 15, ylab = NULL) {
  # Build newdata grid
  newdata <- tibble(mean_temp = temps)
  time_grid <- seq(0, days, by = 1)
  
  # flexsurv::summary returns a list per newdata row
  pred_list <- summary(model, newdata = newdata, type = "survival", t = time_grid, ci = TRUE)
  
  
  pred_df <- map2_dfr(pred_list, newdata$mean_temp, ~ tibble(
    time = .x$time,
    surv = .x$est,
    lcl = .x$lcl,
    ucl = .x$ucl,
    mean_temp = .y
  ))
  ggplot(pred_df, aes(x = time, y = surv, color = factor(mean_temp))) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = PLOT_PALETTE, name = "Mean Temperature", labels = paste0(temps, " °C")) +
    labs(x = "Days after entering stage", y = ylab, title = stage_title) +
    theme_bw(base_size = 14)
}


# Fit AFT models and generate curves for all stages


fit_single_aft <- function(stage_name) {
  dat <- subset_data(stage_name)
  if (nrow(dat) < 5) return(NULL)
  flexsurvreg(Surv(duration, status) ~ mean_temp, data = dat, dist = "weibull")
}


m_first <- fit_single_aft("first")
m_second <- fit_single_aft("second")
m_third <- fit_single_aft("third")
m_fourth <- fit_single_aft("fourth")
m_fifth <- fit_single_aft("fifth")
m_pupa <- fit_single_aft("pupa")


p1 <- plot_aft_survival_curves(m_first, "First instar", ylab = "Probability of remaining in stage")
p2 <- plot_aft_survival_curves(model_second_instar, "Second instar")
p3 <- plot_aft_survival_curves(model_third_instar, "Third instar")
p4 <- plot_aft_survival_curves(model_fourth_instar, "Fourth instar", ylab = "Probability of remaining in stage")
p5 <- plot_aft_survival_curves(model_fifth_instar, "Fifth instar")

legend <- get_legend(p1)

cowplot::plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"), 
  p3 + theme(legend.position = "none"), 
  p4 + theme(legend.position = "none"),
  p5 + theme(legend.position = "none"),
  legend)
