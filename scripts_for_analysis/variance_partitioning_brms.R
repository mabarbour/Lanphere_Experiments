## NOTE THAT PREDICTORS SHOULD BE CENTERED PRIOR TO ANALYSIS.

get_variance_wind_brms <- function(brms_model, data, RE_rows) {
  require(broom)
  tidy.brm <- tidy(brms_model) %>% select(-std.error)
  
  RE.data <- tidy.brm[RE_rows, ]
  
  int <- filter(tidy.brm, term == "b_Intercept")$estimate
  bWind <- filter(tidy.brm, term == "b_Wind.Exposure1")$estimate
  #up95.Wind <- filter(tidy.brm, term == "b_Wind.Exposure1")$upper
  #low95.Wind <- filter(tidy.brm, term == "b_Wind.Exposure1")$lower
  wind.fixef.design <- tidy(model.matrix(~Wind.Exposure, data = data))
  
  avg.wind <- sd(int*wind.fixef.design$X.Intercept. + bWind*wind.fixef.design$Wind.Exposure1)
  #upper.wind <- sd(int*wind.fixef.design$X.Intercept. + up95.Wind*wind.fixef.design$Wind.Exposure1)
  #lower.wind <- sd(int*wind.fixef.design$X.Intercept. + low95.Wind*wind.fixef.design$Wind.Exposure1)
  
  FE.data <- data.frame(term="sd_Wind.Exposure", estimate=avg.wind, lower=NA, upper=NA) #lower = lower.wind, upper = upper.wind)
  
  variance.data <- bind_rows(FE.data, RE.data) #%>%
    #mutate(percent_variance = estimate^2/(sum(estimate^2))) %>%
    #select(term, percent_variance, estimate, lower, upper)
  
  return(variance.data)
}

