
# estimate a simple model to predict the number of confirmed cases a week ahead. 
# The model is a random coefficients exponential growth model. To improve forecasts,
# I weight recent observations more heavily, with the decay rate chosen to minimize
# 7-day-ahead forecast error. 

setwd("~/Dropbox/random/covid19/")
options(stringsAsFactors = FALSE)


set.seed(5103826)


library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)
library(lme4)
library(assertthat)



# Download and clean data -------------------------------------------------


# data from github repo CSSEGISandData/COVID-19
case_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv"
dat = read.csv(case_url)
names(dat)[grepl("^X", names(dat))] = gsub("X", "date_", names(dat)[grepl("^X", names(dat))])

# reshape data to be in long format
dat = dat %>% 
  rename(state = Province.State,
         country = Country.Region,
         lat = Lat, long = Long)

dat = dat %>% 
  select(country, starts_with("date_")) %>% 
  group_by(country) %>% 
  summarise_all(sum) %>% 
  pivot_longer(-country, names_to = "date", values_to = "cases")

# format date
dat = dat %>% 
  mutate(date = gsub("^date_", "", date)) %>% 
  mutate(date = gsub(".20$", ".2020", date,fixed = TRUE)) %>% 
  mutate(date = gsub(".", "-", date, fixed = TRUE)) %>% 
  mutate(date = mdy(date))

# get latest day for which there's data
today = max(dat$date)

# generate log(# cases)
dat = dat %>% 
  mutate(log_cases = ifelse(cases == 0, NA, log(cases))) 

# get number of total cases and days since first case
dat = dat %>% 
  group_by(country) %>% 
  filter(cases > 0) %>% 
  mutate(first_case = min(date,na.rm=TRUE)) %>% 
  mutate(days_since_first_case = date - first_case) %>% 
  mutate(total_cases = max(cases))

# rename some countries
dat = dat %>% 
  ungroup %>% 
  filter(!(as.character(country) %in% c("Others", "Cruise Ship"))) %>% 
  mutate(country = case_when(
    country == "US" ~ "United States",
    country == "UK" ~ "United Kingdom",
    country == "Republic of Ireland" ~ "Ireland", 
    country == "United Arab Emirates" ~ "U.A.E.",
    country == "Iran (Islamic Republic of)" ~ "Iran",
    country == "Republic of Korea" ~ "South Korea",
    country == "Hong Kong SAR" ~ "Hong Kong",
    country == "Taipei and environs" ~ "Taiwan",
    country == "Russian Federation" ~ "Russia",
    country == "Mainland China" ~ "China",
    country == "Korea, South" ~ "South Korea",
    TRUE ~ country
  )) %>% 
  arrange(-total_cases) %>% 
  mutate(country = factor(country, unique(country)))

# only keep countries where the first reported case was at least 10 days ago
dat = dat %>% 
  filter(first_case <= today - 10)


# plot growth in all countries 
ggplot(dat) + 
  aes(x = days_since_first_case, y = log_cases, group = country) + 
  geom_smooth(method = "lm", se = F)




# Tune decay rate ---------------------------------------------------------


# fit model of the form: log cases_it = a_i + b_i Days Since First Case. Fit HLM
# w/ random slopes and intercepts by country. Weight recent observations
# more heavily since they're more useful in determining current trend.
# Observations get exponentially decaying weights with decay rate delta s.t.
# weight for an observation k days ago is delta^k. In order to select delta,
# we'll pretend it's a week ago and try to get the best 7-day-ahead prediction
# for the most recent data we have. I.e., only use training data from 7 or
# more days ago, then find the delta that minimizes the MSE for the prediction
# today.

l = 7
decay = data.frame(delta = seq(.25, 1, .01), 
                   rmse_log = NA_real_, 
                   rmse_level = NA_real_,
                   rmse_pct = NA_real_,
                   mae_level = NA_real_,
                   mae_log = NA_real_,
                   mae_pct = NA_real_) # delta candidates
median_today = dat %>% distinct(country, total_cases) %>% .$total_cases %>% median

# training data is data from no more than 7 days ago. Also discard
# countries that didn't have at least 10 days' worth of data.
traindat = dat %>%
  filter(date <= max(date) - l) %>%
  mutate(days_ago = as.numeric(max(date) - date)) %>%
  filter(total_cases >= median_today & first_case <= max(date) - l - 10)

# hold out data is just the data from today
testdat = dat %>%
  filter(date == max(date) & country %in% traindat$country)

for (d in decay$delta){
  mod = lmer(formula = log_cases ~ days_since_first_case + (days_since_first_case | country),
             data = traindat, control = lmerControl(optimizer ="Nelder_Mead"),
             weights = d^traindat$days_ago)
  preds = predict(mod, newdata = testdat)
  
  err_level = testdat$total_cases - exp(preds)
  err_log   = testdat$log_cases - preds
  err_pct   = 100 * err_level / testdat$total_cases
  
  rmse_log = sqrt(mean(err_log^2))
  rmse_level = sqrt(mean(err_level^2))
  rmse_pct   = sqrt(mean(err_pct^2))
  mae_level = mean(abs(err_level))
  mae_log = mean(abs(err_log))
  mae_pct = mean(abs(err_pct))
  
  decay$rmse_log[decay$delta == d] = rmse_log
  decay$rmse_level[decay$delta == d] = rmse_level
  decay$rmse_pct[decay$delta == d] = rmse_pct
  decay$mae_level[decay$delta == d] = mae_level
  decay$mae_log[decay$delta == d] = mae_log
  decay$mae_pct[decay$delta == d] = mae_pct
}

decay_long = decay %>% 
  pivot_longer(-delta, names_to = "stat") %>% 
  filter(grepl("pct", stat))


ggplot(subset(decay_long, delta >=.5)) + 
  aes(x = delta, y = value, colour = stat) + 
  geom_line()



# plot MSE vs. exponential decay weights for holdout data
delta.opt = decay$delta[which.min(decay$rmse_pct)]
assert_that(abs(delta.opt - .66) < 1e4, msg = "the delta.opt is hard coded in the plot below. update it.")

decay_long %>% 
  filter(grepl("pct", stat), 
         delta >= .58) %>% 
ggplot() + 
  aes(y = value, x = delta, lty = stat) +
  geom_line() + 
  annotate(geom = "text", x = .87, y = 100, label = "Root mean squared\nprediction error", hjust = 0) +
  annotate(geom = "text", x = .8, y = 60, label = "Mean absolute\nprediction error", hjust = 0) +
  geom_vline(xintercept = delta.opt, lty = 2) + 
  annotate(geom = "text", x = delta.opt+.01, y = 100, label = expression(delta^opt==0.66), hjust = 0) + 
  guides(lty = FALSE) + 
  labs(x = expression(paste("Decay rate, ", delta)), y = "Prediction Error as Percentage of True Value") + 
  theme_minimal() + 
  scale_x_continuous(breaks = seq(0, 1, .1)) + 
  ggtitle("Prediction Error by Choice of Decay Rate")

fname = sprintf("figs/decay_rate_%s_%s.png", month.name[month(today)], day(today))
ggsave(fname, width = 6, height=4, units = "in")
fname2 = "figs/decay_rate_current.png"
ggsave(fname2, width = 6, height=4, units = "in")

# save rmse data
decay %>% saveRDS(file = "out/decay_rate.rds")



# Estimate model ----------------------------------------------------------


dat$days_ago = as.numeric(today - dat$date)
mod = lmer(formula = log_cases ~ days_since_first_case + (days_since_first_case | country), 
           data = dat, control = lmerControl(optimizer ="Nelder_Mead"),
           weights = delta.opt^dat$days_ago)



# Make predictions --------------------------------------------------------

# predict today + the next 10 days
newdat = expand.grid(country = unique(dat$country),
                     date = seq(today, today + 7, by = "1 day")) 

preddat = bind_rows(dat, newdat) %>% 
  group_by(country) %>% 
  mutate(first_case = min(first_case, na.rm=TRUE), 
         total_cases = max(total_cases, na.rm=TRUE)) %>% 
  mutate(days_since_first_case = date - first_case) %>% 
  mutate(future = date > today)
  

preddat$prediction_log = predict(mod, newdata = preddat)
preddat$prediction_log = ifelse(preddat$prediction_log <= 0, NA, preddat$prediction_log)

ggplot(subset(preddat, total_cases >= 50)) + 
  aes(x = date) + 
  geom_line(aes(y = log_cases, colour = "Observed")) + 
  geom_line(aes(y = prediction_log, colour = "Model fit")) + 
  facet_wrap(~country)




# Prepare data for plotting -----------------------------------------------

# get the model coefficients into a df for labeling
mod.coefs = coef(mod)$country
mod.coefs$country = rownames(mod.coefs)
mod.coefs$country = factor(mod.coefs$country, levels(dat$country))
mod.coefs = mod.coefs %>% 
  rename(intercept = `(Intercept)`,
         slope = days_since_first_case)

# bind together coefficients and data
preddat = left_join(preddat, mod.coefs, by = "country")

# get each country's forecast for a week from now
forecast = preddat %>% 
  filter(date == today + 7) %>% 
  select(country, forecast_7_day = prediction_log) %>% 
  mutate(forecast_7_day = exp(forecast_7_day))
  
preddat = left_join(preddat, forecast, by = "country")
preddat$label = with(preddat, paste0(
  country, 
  sprintf(
    "\n(%s days; %sk cases)", 
      round( log(2) / slope, 1), 
      round(forecast_7_day/1e3, 1)
  )
))
preddat$label = factor(preddat$label, unique(preddat$label))


# Only plotting a subset of countries - keep it at 29 so it's 5x6 with one
# panel for a legend.
countries_to_keep = unique(subset(preddat, country != "" & total_cases >= 50)$country)
if (length(countries_to_keep) > 29) {
  countries_to_keep = countries_to_keep[1:29]
}
plotdf = subset(preddat, country %in% countries_to_keep)



# Plot results ------------------------------------------------------------

pl = ggplot(plotdf) + 
  aes(x = date, y = prediction_log) + 
  geom_vline(xintercept = today, colour = "grey70", lty = 2) +
  geom_line(aes(lty = future, colour = "Model trend")) + 
  geom_line(aes(y = log_cases, colour = "Confirmed cases")) +
  facet_wrap(~label, drop = TRUE) + 
  scale_x_date(date_labels = "%b %e", breaks = seq(ymd("2020-01-01"), today + 7, by = "15 days")) + 
  scale_y_continuous(breaks = c(log(10), log(1e3), log(1e5)), 
                     labels = c("10", "1k", "100k")) + 
  coord_cartesian(xlim = c(today - 14, today + 7), ylim = c(0, log(1e6))) + 
  theme_minimal() + 
  guides(lty = FALSE) + 
  scale_colour_manual(name = NULL, values = c("Confirmed cases" = scales::muted("red"),
                                              "Model trend" = scales::muted("blue"))) + 
  theme(panel.spacing = unit(1.5, "mm"), 
        plot.caption = element_text(hjust = 0, size = 8), 
        legend.position = c(.93, .08), 
        plot.margin = margin(10, 20, 10, 10),
        strip.text = element_text(size = 7))  + 
  labs(x = NULL, y = "Cumulative number of cases (log scale)",
       title = paste0("Confirmed COVID-19 Cases as of ", month.name[month(today)], " ", day(today)), 
       subtitle = paste("Estimated doubling time and projected number of cases on", month.name[month(today + 7)], day(today+7),
                        "listed in parentheses."))

cap = paste("Model is a random coefficients exponential growth model that takes the form: log(# Cases) = a + b x Days Since First Case,",
            "where a and b\nvary by country.",
            "Doubling time is given by log(2)/b.",
            "Observations from k days ago are given weight d^k with d =", delta.opt, "chosen to minimize the\n7-day-ahead forecast error.",
            "Data source: GitHub.com/CSSEGISandData/COVID-19. Graphic by Will Marble (@wpmarble).")

# save plot
fname = paste0("figs/cases_", month.name[month(today)], "_", day(today))
ggsave(pl+labs(caption = cap), filename = paste0(fname,".pdf"), width=8,height=6)
ggsave(pl+labs(caption = cap), filename = paste0(fname,".png"), width=8,height=6, units = "in")  


fname2 = "figs/cases_current"
ggsave(pl+labs(caption = cap), filename = paste0(fname2,".pdf"), width=8,height=6)
ggsave(pl+labs(caption = cap), filename = paste0(fname2,".png"), width=8,height=6, units = "in") 

fname3 = "figs/cases_current_nocaption.png"
ggsave(pl, filename =fname3, width=8,height=6, units = "in") 




# Save projections --------------------------------------------------------


saveRDS(preddat, file = "out/clean_data_projections.rds")
saveRDS(mod, file = "out/model_results.rds")
