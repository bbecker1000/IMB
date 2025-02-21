library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)
library(ggplot2)
library(survival)

dd_data <- read_csv(here::here("data", "dd_data.csv")) %>% 
  mutate(
    imb_id = as.factor(imb_id),
    host_plant = as.factor(host_plant),
    collection_stage = as.factor(collection_stage),
    overall_survival = as.factor(overall_survival),
    stage = as.factor(stage),
    start = ymd(start),
    end = ymd(end)
  ) %>% 
  mutate_at(vars(contains("stage")), 
            factor,
            levels = c("egg", "first", "second", "third", "fourth", "fifth", "pupa"),
            ordered = TRUE)


ggplot(data = dd_data, aes(x = stage)) +
  geom_bar(stat = "count", aes(fill = survival)) +
  theme_bw()

ggplot(data = dd_data, aes(x = collection_stage)) +
  geom_bar(stat = "count", aes(fill = survival)) +
  theme_bw()

ggplot(data = dd_data, aes(x = total_degree_days)) +
  geom_density(aes(fill = stage), alpha = 0.4) +
  theme_bw()

ggplot(data = dd_data, aes(x = duration)) +
  geom_density(aes(fill = stage), alpha = 0.4) +
  theme_bw()

cor(dd_data$duration, dd_data$total_degree_days)

dd_data %>% count(duration)
