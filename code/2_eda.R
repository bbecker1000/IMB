library(dplyr)
library(tidyverse)
library(here)
library(lubridate)
library(forcats)
library(ggplot2)
library(survival)
library(climate)

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

# plot of mean temperatures for each stage
mean_temp_plot <- ggplot(data = temp_data, aes(x = mean_temp, fill = stage)) +
  geom_density(alpha = 0.4) +
  theme_bw()
mean_temp_plot

# plot of mean temperatures for each stage
mean_temp_plot <- ggplot(data = temp_data %>% filter(stage == "pupa"), aes(x = mean_temp, fill = stage2)) +
  geom_density(alpha = 0.4) +
  theme_bw()
mean_temp_plot

# plot of degree days by year
temps <- station_temps_cleaned %>% 
  mutate(year = year(date),
         month = month(date),
         day0 = mdy(paste0("01/01/", year)),
         dayofYear = as.integer(date - day0),
         degree_days = if_else(degree_days == "NaN", 0, degree_days),
         year = as.factor(year)) %>% 
  filter(year != "2024") %>% 
  group_by(year) %>% 
  mutate(cum_dd = cumsum(degree_days))

ggplot(data = temps, aes(x = dayofYear, y = cum_dd)) +
      geom_line(color = "darkgreen", linewidth = 1.3) +
  facet_wrap(~year) +
  theme_bw() +
  labs(x = "day of year (0 = Jan. 1st)", y = "cumulative degree days")


dd_data %>% group_by(imb_id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(overall_survival)


ggplot(data = data, aes(x = stage)) +
  geom_bar(stat = "count", aes(fill = survival)) +
  theme_bw()

ggplot(data = data %>% filter(stage != "pupa"), aes(x = total_degree_days)) +
  geom_density(aes(fill = stage), alpha = 0.4) +
  theme_bw()

ggplot(data = data %>% filter(stage != "pupa"), aes(x = duration)) +
  geom_density(aes(fill = stage), alpha = 0.4) +
  theme_bw()

cor(dd_data$duration, dd_data$total_degree_days)

ggplot(data = data %>% filter(stage != "pupa"), aes(x = duration, y = total_degree_days)) +
  geom_point() +
  geom_smooth(method = "glm")

ggplot(data = data %>% filter(stage == "pupa"), aes(x = duration, y = total_degree_days)) +
  geom_point() +
  geom_smooth(method = "glm")



# reading in temperature data from rearing rooms. deleting from data prep because we're using station temps
threshold <- 10 # really not sure what the temperature threshold is supposed to be...
# cannot find minimum development threshold for IMB, but i found a source that generally says the min threshold for lots
# of butterflies is 15

temp_data <- read_csv(here::here("data", "continuous_temps.csv")) %>% 
  mutate(location = `...1`,
         date = ymd(date),
         time = hms(time),
         temp = (temp.f - 32) * (5/9)) %>% 
  filter(!str_detect(location, "mystery")) %>% 
  select(date, time, temp) %>% 
  group_by(date) %>% 
  summarise(min = min(temp),
            max = max(temp)) %>% 
  ungroup() %>% 
  mutate(degree_days = pmax(0, ((max+min)/2) - threshold))


# testing correlations between rearing room and weather station data to see if extrapolation is reasonable
# requires getting weatehr station temps by running data prep first

temp_comparisons <- left_join(temp_data, station_temps_cleaned, 
                              by = join_by(date == date),
                              suffix = c("_rearing_room", "_weather_station"))

corr <- glm(degree_days_weather_station ~ degree_days_rearing_room,
            data = temp_comparisons)
summary(corr)

c <- cor(temp_comparisons$degree_days_rearing_room, temp_comparisons$degree_days_weather_station)
c

ggplot(data = temp_comparisons, aes(x = degree_days_weather_station, y = degree_days_rearing_room)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "glm") +
  theme_bw()


