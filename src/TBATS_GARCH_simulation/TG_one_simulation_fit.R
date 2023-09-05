library(rugarch)
library(forecast)

library(here)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
i_am('src/TBATS_GARCH_simulation/TG_one_simulation_fit.R')

source(here('src','new_forecast','fitTBATS.R'))
source(here('src','new_forecast','tbats.R'))
source(here('src','new_forecast','forecastTBATS.R'))

attach(loadNamespace('forecast'), name = 'forecast_all')

# Load data
simulation_data <- read.csv(file=here('data','derived','TBATS_GARCH_simulation','TG_one_simulation_data.csv'))
observed <- simulation_data$observed
periods <- c(7, 30.42, 100)

# Fit a TBATS model
tbats_model <- fitSpecificTBATS(observed, use.box.cox=FALSE, use.beta=TRUE, use.damping=FALSE, 
                                seasonal.periods=periods, k.vector=c(1,1,1),
                                bc.lower=0, bc.upper=1, biasadj=FALSE)
# Fit a TBATS model with ARMA errors
tbats_arma_model <- fitSpecificTBATS(observed, use.box.cox=FALSE, use.beta=TRUE, use.damping=FALSE, 
                                seasonal.periods=periods, k.vector=c(1,1,1),
                                ar.coefs=numeric(1), ma.coefs=numeric(1), 
                                bc.lower=0, bc.upper=1, biasadj=FALSE)
# Fit a TBATS model with GARCH errors   
tbats_garch_model <- fitSpecificTBATS(observed, use.box.cox=FALSE, use.beta=TRUE, use.damping=FALSE, 
                                    seasonal.periods=periods, k.vector=c(1,1,1), 
                                    v=1, u=1, bc.lower=0, bc.upper=1, biasadj=FALSE)
# Fit a TBATS model with ARMA and GARCH errors
tbats_arma_garch_model <- fitSpecificTBATS(observed, use.box.cox=FALSE, use.beta=TRUE, use.damping=FALSE, 
                                    seasonal.periods=periods, k.vector=c(1,1,1), 
                                    ar.coefs=numeric(1), ma.coefs=numeric(1), 
                                    v=1, u=1, bc.lower=0, bc.upper=1, biasadj=FALSE)    
# Fit a TBATS model first, and then fit a GARCH model
spec1 <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
                    distribution = 'norm')
tbats_garch_model2 <- ugarchfit(spec = spec1, data = tbats_model$errors)                                             
# Fit a TBATS model with ARMA errors first, and then fit a GARCH model
spec2 <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
                    distribution = 'norm')
tbats_arma_garch_model2 <- ugarchfit(spec = spec2, data = tbats_arma_model$errors)
# Fit a TBATS model first, and then fit a ARMA + GARCH model
spec3 <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 1), include.mean=FALSE),
                    distribution = 'norm')
tbats_arma_garch_model3 <- ugarchfit(spec = spec3, data = tbats_model$errors)

# Print the AIC of the above models
tbats_model$AIC
tbats_arma_model$AIC
tbats_garch_model$AIC
tbats_arma_garch_model$AIC
-2*likelihood(tbats_garch_model2) + tbats_garch_model$AIC-tbats_garch_model$lik
-2*likelihood(tbats_arma_garch_model2) + tbats_arma_garch_model$AIC-tbats_arma_garch_model$lik
-2*likelihood(tbats_arma_garch_model3) + tbats_arma_garch_model$AIC-tbats_arma_garch_model$lik

model1 <- tbats.components(tbats_model)
model2 <- tbats.components(tbats_arma_model)
model3 <- tbats.components(tbats_garch_model)
model4 <- tbats.components(tbats_arma_garch_model)
components1 <- list('level'=model1[,2], 'slope'=model1[,3], 'seasonal1'=model1[,4], 'seasonal2'=model1[,5], 'seasonal3'=model1[,6])
components2 <- list('level'=model2[,2], 'slope'=model2[,3], 'seasonal1'=model2[,4], 'seasonal2'=model2[,5], 'seasonal3'=model2[,6])
components3 <- list('level'=model3[,2], 'slope'=model3[,3], 'seasonal1'=model3[,4], 'seasonal2'=model3[,5], 'seasonal3'=model3[,6])
components4 <- list('level'=model4[,2], 'slope'=model4[,3], 'seasonal1'=model4[,4], 'seasonal2'=model4[,5], 'seasonal3'=model4[,6])
write.csv(components1, file=here('data','derived','TBATS_GARCH_simulation','TG_one_simulation_fit_components1.csv'))
write.csv(components2, file=here('data','derived','TBATS_GARCH_simulation','TG_one_simulation_fit_components2.csv'))
write.csv(components3, file=here('data','derived','TBATS_GARCH_simulation','TG_one_simulation_fit_components3.csv'))
write.csv(components4, file=here('data','derived','TBATS_GARCH_simulation','TG_one_simulation_fit_components4.csv'))