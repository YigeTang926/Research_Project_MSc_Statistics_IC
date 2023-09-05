compute_iteration <- function(i) {
    set.seed(i)
    simulate_arma_garch <- function (n, mu, phi, theta, omega, alpha, beta) {
        p <- length(phi)
        q <- length(theta)
        v <- length(alpha)
        u <- length(beta)
        # simulate GARCH process
        z <- rnorm(n)
        e <- rep(0,n)
        s <- rep(0,n)
        for (i in (max(u,v)+1):n) {
            s[i] <- sqrt(omega + sum(alpha * rev(e[(i-v):(i-1)]^2)) + sum(beta * rev(s[(i-u):(i-1)]^2)))
            e[i] <- s[i] * z[i]
        }
        # simulate ARMA process
        x <- rep(0,n)
        for (i in (max(p,q)+1):n) {
            x[i] <- mu + sum(phi * rev(x[(i-p):(i-1)])) + sum(theta * rev(e[(i-q):(i-1)]))
        }
        return(x)
    }
    n  <- 1000
    # Trend
    trend <- 0.05*(1:n)
    # Seasonal
    periods = c(7,30.42,100)
    amptitudes = c(3,3,3)
    seasonal_1 <- amptitudes[1]*sin(2*pi*(1:n)/periods[1])
    seasonal_2 <- amptitudes[2]*cos(2*pi*(1:n)/periods[2])
    seasonal_3 <- amptitudes[3]*cos(2*pi*(1:n)/periods[3] + pi/4)
    seasonal <- seasonal_1 + seasonal_2 + seasonal_3
    # Simulate ARMA(1,1)+GARCH(1,1) process
    mu = 0
    phi <- c(0.8)
    theta <- c(0.5)
    omega <- 1
    alpha <- c(0.5)
    beta <- c(0.4)
    error <- simulate_arma_garch(n = n, mu = mu, phi = phi, theta = theta, omega = omega, alpha = alpha, beta = beta)
    observed <- trend + seasonal + error 

    tryCatch({
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
    AIC1 <- -2*likelihood(tbats_garch_model2) + tbats_garch_model$AIC-tbats_garch_model$lik                                        
    # Fit a TBATS model with ARMA errors first, and then fit a GARCH model
    spec2 <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)),
                        mean.model = list(armaOrder = c(0, 0), include.mean=FALSE),
                        distribution = 'norm')
    tbats_arma_garch_model2 <- ugarchfit(spec = spec2, data = tbats_arma_model$errors)
    AIC2 <- -2*likelihood(tbats_arma_garch_model2) + tbats_arma_garch_model$AIC-tbats_arma_garch_model$lik
    # Fit a TBATS model first, and then fit a ARMA + GARCH model
    spec3 <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)),
                        mean.model = list(armaOrder = c(1, 1), include.mean=FALSE),
                        distribution = 'norm')
    tbats_arma_garch_model3 <- ugarchfit(spec = spec3, data = tbats_model$errors)
    AIC3 <- -2*likelihood(tbats_arma_garch_model3) + tbats_arma_garch_model$AIC-tbats_arma_garch_model$lik
    
    # Store results in the results matrix
    aic <- c(tbats_model$AIC, tbats_arma_model$AIC, tbats_garch_model$AIC, tbats_arma_garch_model$AIC, AIC1, AIC2, AIC3)
    parameters1 <- c(tbats_arma_garch_model$ar.coefficients, tbats_arma_garch_model$ma.coefficients, 
                    tbats_arma_garch_model$mu, tbats_arma_garch_model$arch.coefficients, tbats_arma_garch_model$garch.coefficients)
    parameters2 <- c(tbats_arma_model$ar.coefficients, tbats_arma_model$ma.coefficients,
                    coef(tbats_arma_garch_model2)[1], coef(tbats_arma_garch_model2)[2], coef(tbats_arma_garch_model2)[3])
    parameters3 <- c(coef(tbats_arma_garch_model3)[1], coef(tbats_arma_garch_model3)[2], 
                    coef(tbats_arma_garch_model3)[3], coef(tbats_arma_garch_model3)[4], coef(tbats_arma_garch_model3)[5])
    c(aic, parameters1, parameters2, parameters3)
    }, error = function(e) {
        numeric(22)
    })
}

# Parallel Calculating
library(foreach)
library(doParallel)
num_cores <- detectCores()
cl <- makeCluster(num_cores - 4)
registerDoParallel(cl)

stopImplicitCluster()

library(doSNOW)
registerDoSNOW(cl)

# progress bar
iterations <- 1000
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)

opts <- list(progress = progress)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
result <- foreach(j = 1:iterations, .combine = rbind, .options.snow = opts, .packages=c('forecast', 'rugarch', 'here')) %dopar%
{   
    i_am('src/TBATS_GARCH_simulation/TG_1000_simulations_parallel.R')
    source(here('src','new_forecast','fitTBATS.R'))
    source(here('src','new_forecast','tbats.R'))
    source(here('src','new_forecast','forecastTBATS.R'))
    attach(loadNamespace('forecast'), name = 'forecast_all')
    compute_iteration(j)
}

close(pb)
stopCluster(cl)

# Delete the rows with all zero values and unreasonable values
all_zero_rows <- apply(result, 1, function(row) all(row == 0))
results <- result[!all_zero_rows, ]
all_unreasonable_rows <- apply(results[,8:22], 1, function(row) any(row > 10))
results <- results[!all_unreasonable_rows, ]

# Calculate the means of the results
mean <- round(colMeans(results),4)
aic_mean <- mean[1:7]
parameters1_mean <- mean[8:12]
parameters2_mean <- mean[13:17]
parameters3_mean <- mean[18:22]
# Calculate the standard deviations of the results
std <- round(apply(results, 2, sd),4)
aic_std <- std[1:7]
parameters1_std <- std[8:12]
parameters2_std <- std[13:17]
parameters3_std <- std[18:22]

# Save the results
write.csv(results, file = here('data', 'derived', 'TBATS_GARCH_simulation', 'TG_1000_simulations_results.csv'), row.names = FALSE)
write.csv(mean, file = here('data', 'derived', 'TBATS_GARCH_simulation', 'TG_1000_simulations_mean.csv'), row.names = FALSE)
write.csv(std, file = here('data', 'derived', 'TBATS_GARCH_simulation', 'TG_1000_simulations_std.csv'), row.names = FALSE)