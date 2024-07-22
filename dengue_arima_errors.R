

library(readxl)
library(tidyverse)
library(fpp3)
library(tsibble)
library(ggplot2)



base_dengue <- read_excel("dengue_dataset.xlsx")

### tsibble
base_dengue <- base_dengue  %>%
  mutate(Month = as.Date(
    paste(ano, mes, "01", sep = "-"),
    format = "%Y-%m-%d")) %>% 
  select(-ano, -mes) %>% 
  mutate(Month = yearmonth(Month)) %>% 
  as_tsibble(index= Month) 

# removing pre-2013 data and creating dummy

base_dengue_pos2013 <- base_dengue[181:300,]

base_dengue_pos2013_correc <- base_dengue_pos2013 %>% mutate(D1000 = 0)

base_dengue_pos2013_correc <- base_dengue_pos2013_correc %>%
  mutate(D1000 = ifelse(casosmais1_por_100k >= 1000, 1, D1000))

rm(base_dengue, base_dengue_pos2013)



# acf

# dengue cases
base_dengue_pos2013_correc |>
  gg_tsdisplay(log(casosmais1_por_100k),
               plot_type='partial', lag=24) +
  labs(title="Log (dengue cases plus one per 100,000 population)", y="")

# temperature
base_dengue_pos2013_correc |>
  gg_tsdisplay(log(tmed),
               plot_type='partial', lag=24) +
  labs(title="Log (temperature)", y="")

# precipitation
base_dengue_pos2013_correc |>
  gg_tsdisplay(log(chuva),
               plot_type='partial', lag=24) +
  labs(title="Log (precipitation)", y="")



# plotting logarithm of series

base_graf <- base_dengue_pos2013_correc

base_graf <- base_graf %>% rename('Cases (+1) per 100,000 pop.' = casosmais1_por_100k,
                                  'Mean temperature' = tmed,
                                  'Precipitation' = chuva)



base_graf %>%
  pivot_longer(c(`Cases (+1) per 100,000 pop.`, `Mean temperature`, Precipitation), names_to = "var", values_to = "value") %>%
  mutate(value = ifelse(value > 0, log(value), NA)) %>%
  ggplot(aes(x = Month, y = value)) +
  geom_line() +
  facet_grid(vars(var), scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "",
       y = "Log (value)", x = "Month")




# model 1

# pure ARIMA

reg_log1 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA(log(casosmais1_por_100k),
                      greedy = FALSE,
                      stepwise = FALSE,
                      approximation = FALSE
  )) %>% 
  report(reg_log1) 

reg_log1 |> gg_tsresiduals()

augment(reg_log1) |> features(.innov, ljung_box, dof=4, lag=16)
augment(reg_log1) |> features(.innov, ljung_box, dof=4, lag=24)



resid_model1 <- reg_log1 %>%
  residuals() %>%
  as_tibble()


resid_model1_ts <- ts(resid_model1[3], frequency = 12, start = c(2013,1))

# heteroskedasticity

FinTS::ArchTest(resid_model1_ts)


# normality

# Shapiro-Wilk 

shapiro.test(resid_model1_ts)




# model 2

# Model 2: Controlling for contemporaneous precipitation and temperature

reg_log2 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA(log(casosmais1_por_100k) ~ 
                        log(tmed) + log(chuva),
                      greedy = FALSE,
                      stepwise = FALSE
  )) %>% 
  report(reg_log2) 

reg_log2 |> gg_tsresiduals()

augment(reg_log2) |> features(.innov, ljung_box, dof=4, lag=16)
augment(reg_log2) |> features(.innov, ljung_box, dof=4, lag=24)



resid_model2 <- reg_log2 %>%
  residuals() %>%
  as_tibble()


resid_model2_ts <- ts(resid_model2[3], frequency = 12, start = c(2013,1))

# heteroskedasticity

FinTS::ArchTest(resid_model2_ts)


# normality

# Shapiro-Wilk 

shapiro.test(resid_model2_ts)




# model 3

# Model 3: Controlling for contemporaneous precipitation and temperature, plus a one-period lag

reg_log3 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA( log(casosmais1_por_100k) ~ 
                         log(tmed) + log(lag(tmed, 1))  +
                         log(chuva) + log(lag(chuva, 1)),
                       greedy = FALSE,
                       stepwise = FALSE,
                       approximation = FALSE
  )) %>% 
  report(reg_log3) 

reg_log3 |> gg_tsresiduals()

augment(reg_log3) |> features(.innov, ljung_box, lag=18, dof = 6)
augment(reg_log3) |> features(.innov, ljung_box, lag=24, dof = 6)




resid_model3 <- reg_log3 %>%
  residuals() %>%
  as_tibble()


resid_model3_ts <- ts(resid_model3[3], frequency = 12, start = c(2013,1))

# heteroskedasticity

FinTS::ArchTest(resid_model3_ts)


# normality

# Shapiro-Wilk 

shapiro.test(resid_model3_ts)




# model 4

# Model 4: Controlling for contemporaneous precipitation and temperature, plus a two-period lag

reg_log4 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA( log(casosmais1_por_100k) ~ 
                         log(tmed) + log(lag(tmed, 1)) + log(lag(tmed, 2)) +
                         log(chuva) + log(lag(chuva, 1)) + log(lag(chuva, 2)),
                       greedy = FALSE,
                       stepwise = FALSE,
                       approximation = FALSE
  )) %>% 
  report(reg_log4) 

reg_log4 |> gg_tsresiduals()

# report(reg_log1)

augment(reg_log4) |> features(.innov, ljung_box, lag=18, dof = 3)
augment(reg_log4) |> features(.innov, ljung_box, lag=24, dof = 3)



# ts object

resid_model4 <- reg_log4 %>%
  residuals() %>%
  as_tibble()

resid_model4_ts <- ts(resid_model4[3], frequency = 12, start = c(2013,1))

# heteroskedasticity

FinTS::ArchTest(resid_model4_ts)


# normality

# Shapiro-Wilk 

shapiro.test(resid_model4_ts)





# model 5 (final model)

# Model 5 (final model): Controlling for contemporaneous precipitation and temperature, plus a one-period lag, plus a dummy for the three months which have >1000 cases per 100,000 population

reg_log5 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA( log(casosmais1_por_100k) ~ 
                         log(tmed) + log(lag(tmed, 1)) +
                         log(chuva) + log(lag(chuva, 1)) +
                         D1000,
                       greedy = FALSE,
                       stepwise = FALSE,
                       approximation = FALSE
  )) %>% 
  report(reg_log5) 

reg_log5 |> gg_tsresiduals()

# report(reg_log1)

augment(reg_log5) |> features(.innov, ljung_box, lag=18, dof = 6)
augment(reg_log5) |> features(.innov, ljung_box, lag=24, dof = 6)
accuracy(reg_log5, measures = point_accuracy_measures)



# diagnostic tests

# ts object

resid_model5 <- reg_log5 %>%
  residuals() %>%
  as_tibble()

resid_model5_ts <- ts(resid_model5[3], frequency = 12, start = c(2013,1))

# heteroskedasticity

FinTS::ArchTest(resid_model5_ts)


# normality

# Shapiro-Wilk 

shapiro.test(resid_model5_ts)



# manually calculating the SRMSE

residuos_modelo5_sem_NAs <- augment(reg_log5) %>% na.exclude()

desvpad5 <- sd(residuos_modelo5_sem_NAs$casosmais1_por_100k)

rmse5 <- sqrt(mean((residuos_modelo5_sem_NAs$casosmais1_por_100k - residuos_modelo5_sem_NAs$.fitted)^2))

print(rmse5)
# check it with the RMSE provided by 'accuracy' function above - OK

srmse5 <- rmse5/desvpad5

print(srmse5)




# model 6 (final model)

# Model 6 (final model): Controlling for contemporaneous precipitation and temperature, plus a two-period lag, plus a dummy for the three months which have >1000 cases per 100,000 population

reg_log6 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA( log(casosmais1_por_100k) ~ 
                         log(tmed) + log(lag(tmed, 1)) + log(lag(tmed, 2)) +
                         log(chuva) + log(lag(chuva, 1)) + log(lag(chuva, 2)) +
                         D1000,
                       greedy = FALSE,
                       stepwise = FALSE,
                       approximation = FALSE
  )) %>% 
  report(reg_log6) 

reg_log6 |> gg_tsresiduals()

# report(reg_log1)

augment(reg_log6) |> features(.innov, ljung_box, lag=18, dof = 6)
augment(reg_log6) |> features(.innov, ljung_box, lag=24, dof = 6)
accuracy(reg_log6, measures = point_accuracy_measures)



# diagnostic tests

# ts object

resid_model6 <- reg_log6 %>%
  residuals() %>%
  as_tibble()

resid_model6_ts <- ts(resid_model6[3], frequency = 12, start = c(2013,1))

# heteroskedasticity

FinTS::ArchTest(resid_model6_ts)


# normality

# Shapiro-Wilk 

shapiro.test(resid_model6_ts)



# SRMSE

residuos_modelo6_sem_NAs <- augment(reg_log6) %>% na.exclude()

desvpad6 <- sd(residuos_modelo6_sem_NAs$casosmais1_por_100k)

rmse6 <- sqrt(mean((residuos_modelo6_sem_NAs$casosmais1_por_100k - residuos_modelo6_sem_NAs$.fitted)^2))

print(rmse6)
# check it with the RMSE provided by 'accuracy' function above - OK

srmse6 <- rmse6/desvpad6

print(srmse6) 







# appendix: additional models


# model with 3 lags
mod_3lags = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA( log(casosmais1_por_100k) ~ 
                         log(tmed) + log(lag(tmed, 1)) + log(lag(tmed, 2)) + log(lag(tmed, 3)) +
                         log(chuva) + log(lag(chuva, 1)) + log(lag(chuva, 2)) + log(lag(chuva, 3)) +
                         D1000,
                       greedy = FALSE,
                       stepwise = FALSE,
                       approximation = FALSE
  )) %>% 
  report(mod_3lags) 

mod_3lags |> gg_tsresiduals()

# report(reg_log1)

augment(mod_3lags) |> features(.innov, ljung_box, lag=18, dof = 6)
augment(mod_3lags) |> features(.innov, ljung_box, lag=24, dof = 6)



# models with temperature only

mod_temp1 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA( log(casosmais1_por_100k) ~ 
                         log(tmed) + log(lag(tmed, 1)) +
                         D1000,
                       greedy = FALSE,
                       stepwise = FALSE,
                       approximation = FALSE
  )) %>% 
  report(mod_temp1) 

mod_temp1 |> gg_tsresiduals()

augment(mod_temp1) |> features(.innov, ljung_box, lag=18, dof = 3)
augment(mod_temp1) |> features(.innov, ljung_box, lag=24, dof = 3)



mod_temp2 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA( log(casosmais1_por_100k) ~ 
                         log(tmed) + log(lag(tmed, 1)) + log(lag(tmed, 2)) +
                         D1000,
                       greedy = FALSE,
                       stepwise = FALSE,
                       approximation = FALSE
  )) %>% 
  report(mod_temp2) 

mod_temp2 |> gg_tsresiduals()

augment(mod_temp2) |> features(.innov, ljung_box, lag=18, dof = 3)
augment(mod_temp2) |> features(.innov, ljung_box, lag=24, dof = 3)



# models with precipitation only

mod_precip1 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA( log(casosmais1_por_100k) ~ 
                         log(chuva) + log(lag(chuva, 1)) + 
                         D1000,
                       greedy = FALSE,
                       stepwise = FALSE,
                       approximation = FALSE
  )) %>% 
  report(mod_precip1) 

mod_precip1 |> gg_tsresiduals()

augment(mod_precip1) |> features(.innov, ljung_box, lag=18, dof = 6)
augment(mod_precip1) |> features(.innov, ljung_box, lag=24, dof = 6)




mod_precip2 = base_dengue_pos2013_correc %>%
  as_tsibble() %>%
  model(arima = ARIMA( log(casosmais1_por_100k) ~ 
                         log(chuva) + log(lag(chuva, 1)) + log(lag(chuva, 2)) +
                         D1000,
                       greedy = FALSE,
                       stepwise = FALSE,
                       approximation = FALSE
  )) %>% 
  report(mod_precip2) 

mod_precip2 |> gg_tsresiduals()


augment(mod_precip2) |> features(.innov, ljung_box, lag=18, dof = 6)
augment(mod_precip2) |> features(.innov, ljung_box, lag=24, dof = 6)


