library(here)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
i_am('src/data_analysis_czech/seasonal_decomposition_TBATS.R')
library(rugarch)
library(forecast)

source(here('src','new_forecast','fitTBATS.R'))
source(here('src','new_forecast','tbats.R'))
source(here('src','new_forecast','forecastTBATS.R'))

attach(loadNamespace('forecast'), name = 'forecast_all')

path <- here('data','derived','czech')

# Load data
wind_production <- read.csv(file=here(path,'wind_time_series.csv') ,row.names=1)
solar_production <- read.csv(file=here(path,'solar_time_series.csv') ,row.names=1)
wind_periods <- read.csv(file=here(path,'wind_periods.csv') ,row.names=1)
solar_periods <- read.csv(file=here(path,'solar_periods.csv') ,row.names=1)

# Function for seasonal decomposition
seasonal_decomposition <- function(data, period){
  ts <- ts(data, start=c(2012,1), frequency=365.5)
  period <- na.omit(period)
  model <- tbats(ts, use.box.cox=FALSE, use.trend=FALSE, use.arma.errors=TRUE, use.garch.errors=TRUE, seasonal.periods=period)
  if (!is.null(model$seasonal.periods)){
    components <- tbats.components(model)
    print(model$lik/-2)
    return(list(
      'components'=components, 
      'remainder'=components[,1]-rowSums(components[,-1]),
      'order'=c(model$p,model$q,model$v,model$u),
      'k.vector'=model$k.vector,
      'logLik'=model$lik/-2,
      'AIC'=model$AIC      
      )
    )
  } else {
    return(NULL)
  }
}

# Seasonal decomposition (Time Consuming)
output_list_wind <- Map(seasonal_decomposition, wind_production, wind_periods)
output_list_solar <- Map(seasonal_decomposition, solar_production, solar_periods)

# Get the outputs of the output_list_wind
wind_components <- lapply(output_list_wind, function(x) x$components)
wind_remainder <- lapply(output_list_wind, function(x) x$remainder)
wind_order <- lapply(output_list_wind, function(x) x$order)
wind_k <- lapply(output_list_wind, function(x) x$k.vector)
wind_logLik <- lapply(output_list_wind, function(x) x$logLik)
wind_AIC <- lapply(output_list_wind, function(x) x$AIC)

# Save as csv
write.csv(wind_components, file=here(path,'seasonal_decomposition','wind_components.csv'))
write.csv(wind_remainder, file=here(path,'seasonal_decomposition','wind_remainder.csv'))
write.csv(wind_order, file=here(path,'seasonal_decomposition','wind_order.csv'))
write.csv(wind_k, file=here(path,'seasonal_decomposition','wind_k.csv'))
write.csv(wind_logLik, file=here(path,'seasonal_decomposition','wind_logLik.csv'))
write.csv(wind_AIC, file=here(path,'seasonal_decomposition','wind_AIC.csv'))

# Get the outputs of the output_list_wind
solar_components <- lapply(output_list_solar, function(x) x$components)
solar_remainder <- lapply(output_list_solar, function(x) x$remainder)
solar_order <- lapply(output_list_solar, function(x) x$order)
solar_k <- lapply(output_list_solar, function(x) x$k.vector)
solar_logLik <- lapply(output_list_solar, function(x) x$logLik)
solar_AIC <- lapply(output_list_solar, function(x) x$AIC)
# Save as csv
write.csv(solar_components, file=here(path,'seasonal_decomposition','solar_components.csv'))
write.csv(solar_remainder, file=here(path,'seasonal_decomposition','solar_remainder.csv'))
write.csv(solar_order, file=here(path,'seasonal_decomposition','solar_order.csv'))
write.csv(solar_k, file=here(path,'seasonal_decomposition','solar_k.csv'))
write.csv(solar_logLik, file=here(path,'seasonal_decomposition','solar_logLik.csv'))
write.csv(solar_AIC, file=here(path,'seasonal_decomposition','solar_AIC.csv'))