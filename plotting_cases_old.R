
setwd("~/Dropbox/random/covid19/")
options(stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)
library(lme4)


today = as.Date(file.info("COVID19_2020_open_line_list - Hubei.csv")$mtime)

#  Read in data - downloaded from:
# https://docs.google.com/spreadsheets/d/1itaohdPiAeniCXNlntNztZ_oRvjh0HsGuJXUJWET008/edit#gid=0
sheet1 = read.csv("COVID19_2020_open_line_list - Hubei.csv")
sheet2 = read.csv("COVID19_2020_open_line_list - outside_Hubei.csv")



# bind data together
dat = bind_rows(mutate_all(sheet1, as.character), mutate_all(sheet2, as.character))
dat$date = dmy(dat$date_confirmation)
dat = bind_rows(dat, data.frame(country = "China", date = seq(ymd("2020-03-01"), ymd("2020-03-15"), by = "day")))

# tidy data
case_by_country_day = dat %>% 
  filter(!is.na(date)) %>% 
  mutate(country = case_when(
    country != "" ~ country,
    country == "" & province == "Taiwan" ~ "Taiwan"
  )) %>% 
  group_by(country, date) %>% 
  summarise(new_cases = n()) %>% 
  arrange(country, date) %>% 
  ungroup() %>% 
  complete(date, country, fill = list(new_cases = 0)) %>% 
  mutate(new_cases = ifelse(date > today, NA, new_cases)) %>% 
  group_by(country) %>% 
  arrange(country, date) %>% 
  mutate(cumul_cases = cumsum(new_cases),
         cumul_cases_log = log(cumul_cases),
         total_cases = max(cumul_cases, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(country) %>% 
  filter(cumul_cases >= 1 | date > today) %>% 
  mutate(first_case = min(date,na.rm=TRUE)) %>% 
  mutate(days_since_first_case = date - first_case) 

# plot growth (log scale) 
ggplot(case_by_country_day) +
  aes(x = date, y = cumul_cases_log, group = country) + 
  geom_line() + 
  guides(group = FALSE) 

ggplot(case_by_country_day) +
  aes(x = date, y = cumul_cases_log) + 
  geom_line() + 
  guides(group = FALSE) +
  facet_wrap(~country)

ggplot(subset(case_by_country_day, country == "United States")) + 
  aes(x = date, y = cumul_cases) +
  geom_line()
  


# fit a model of the form log(cases_it) = a_i + b_i DaysSinceFirstCase_it
# with random coefficients for a (intercept) and b (slope on time). Corresponds
# to exponential growth model of cases_it = exp(a_i) exp(b_i * DaysSinceFirstCase_it).
fitdat = subset(case_by_country_day, date <= today)
mod = lmer(formula = cumul_cases_log ~ days_since_first_case + (days_since_first_case | country), 
           data = fitdat, control = lmerControl(optimizer ="Nelder_Mead"))
case_by_country_day$prediction_log = predict(mod, newdata = case_by_country_day)


# rename UAE for plotting purposes, include an indicator for whether the date
# is in the future or not
case_by_country_day = case_by_country_day %>% 
  ungroup() %>% 
  arrange(-total_cases) %>% 
  mutate(country2 = car::recode(country, "'United Arab Emirates' = 'U.A.E.'")) %>% 
  mutate(country2 = factor(country2, unique(country2))) %>% 
  mutate(future = ifelse(date > today, TRUE, FALSE)) 


# get the model coefficients into a df for labeling
mod.coefs = coef(mod)$country
mod.coefs$country2 = car::recode(rownames(mod.coefs), "'United Arab Emirates' = 'U.A.E.'")
mod.coefs$country2 = factor(mod.coefs$country2, levels(case_by_country_day$country2))
mod.coefs$label = sprintf(
  "b = %s", 
  round(mod.coefs$days_since_first_case, 3) 
)
mod.coefs$date = ymd("2020-02-15")
mod.coefs$prediction_log = ifelse(mod.coefs$country2 == "China", log(100), log(100000))
mod.coefs$future = NA


# Make the plot - only countries w/ >= 12 confirmed cases
countries_to_keep = unique(subset(case_by_country_day, country != "" & total_cases >= 15)$country2)
if (length(countries_to_keep) > 24) {
  countries_to_keep = countries_to_keep[1:24]
}
plotdf = subset(case_by_country_day, country2 %in% countries_to_keep)
mincases = min(plotdf$total_cases)
pl = ggplot(plotdf) + 
  aes(x = date, y = prediction_log, lty = future) + 
  geom_line(aes(y = cumul_cases_log), colour = "grey30") +
  geom_line(colour = "black") + 
  facet_wrap(~country2, drop = TRUE) + 
  geom_text(
    data = subset(mod.coefs, country2 %in% countries_to_keep),
    size = 2,
    inherit.aes = TRUE,
    aes(label = label)
  ) +
  scale_x_date(date_labels = "%b", breaks = seq(ymd("2020-01-01"), ymd("2020-03-01"), by = "month")) + 
  scale_y_continuous(breaks = c(log(10), log(1e3), log(1e5)), 
                     labels = c("10", "1k", "100k")) + 
  scale_linetype_discrete(name=NULL, label = c("Fitted", "Projection")) +
  theme_minimal() + 
  theme(panel.spacing = unit(1.5, "mm"), 
        plot.caption = element_text(hjust = 0, size = 8), 
        legend.position = c(.92, .08))  + 
  labs(x = NULL, y = "Cumulative number of cases (log scale)",
       title = "Confirmed COVID-19 Cases", 
       caption = paste0("Gray line = confirmed cases. Black line = fitted trend from random coefficients exponential growth\nmodel. Model takes the form log(# Cases) = a + b x Days Since First Case. Equivalently,\n# Cases = exp(a) x exp(b x Days Since First Case). Countries with at least ", 
       mincases, " confirmed cases are shown.\nData source: https://tinyurl.com/s6gsq5y. Graphic by Will Marble (@wpmarble).")) 

fname = paste0("figs/cases_", month(today), "_", day(today))
ggsave(pl, filename = paste0(fname,".pdf"), width=6,height=4)
ggsave(pl, filename = paste0(fname,".png"), width=6,height=4, units = "in")  

fname2 = "figs/cases_current"
ggsave(pl, filename = paste0(fname2,".pdf"), width=6,height=4)
ggsave(pl, filename = paste0(fname2,".png"), width=6,height=4, units = "in")  

