library(rugarch)
library(GNAR)

library(here)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
i_am('src/GNAR_GARCH_simulation/GG_one_simulation_data.R')

source(here('src','new_GNAR','GNARfit.R'))
source(here('src','new_GNAR','logLik.GNARfit.R'))

# Number of observeds
n <- 1000

# Simulate GNAR + GARCH process
simulate_gnar_garch <- function (n, global_alpha, beta, mu, a, b) {
    p <- length(global_alpha)
    s <- lengths(beta)
    v <- length(a[[1]])
    u <- length(b[[1]])
    # Simulate GNAR(p,s) process
    nNodes <- length(fiveNet$edges)
    data_gnar <- GNARsim(n=n, net=fiveNet, alphaParams=lapply(global_alpha, function(x) rep(x,nNodes)),betaParams=beta)
    # simulate GARCH process
    data_garch <- matrix(0,nrow=n,ncol=nNodes)
    for (r in 1:nNodes) {
        z <- rnorm(n)
        e <- rep(0,n)
        s <- rep(0,n)
        for (i in (max(u,v)+1):n) {
            s[i] <- sqrt(mu[r] + sum(a[[r]] * rev(e[(i-v):(i-1)]^2)) + sum(b[[r]] * rev(s[(i-u):(i-1)]^2)))
            e[i] <- s[i] * z[i]
        }
        data_garch[,r] <- e
    }
    data <- data_gnar + data_garch
    return(data)
}

# Parameters
global_alpha <- c(0.3, 0.2)
beta <- list(c(0.1, 0.05),c(0.1))
mu <- rep(0.1,5)
a <- list(c(0.5),c(0.4),c(0.3),c(0.3),c(0.4))
b <- list(c(0.4),c(0.3),c(0.2),c(0.3),c(0.4))

set.seed(111)
iteration <- 1000
matrix <- matrix(NA, nrow=iteration, ncol=8+2*(length(global_alpha)+sum(lengths(beta))))
for (i in 1:iteration) {
    data <- simulate_gnar_garch(n=n,global_alpha=global_alpha, beta=beta, mu=mu, a=a, b=b)
    # Fit a GNAR model
    gnar_model <- GNARfit(vts=data, net=fiveNet, alphaOrder=length(global_alpha),betaOrder=lengths(beta))
    # Fit a GNAR + GARCH model
    gnar_garch_model <- GNARfit(vts=data, net=fiveNet, alphaOrder=length(global_alpha),betaOrder=lengths(beta),archOrder=length(a[1]),garchOrder=length(b[1]))
    # log-likelihood
    logLiks1 <- c(logLik(gnar_model)$ll1, logLik(gnar_model)$ll2, logLik(gnar_model)$ll3)
    logLiks2 <- c(logLik(gnar_garch_model)$ll1, logLik(gnar_garch_model)$ll2, logLik(gnar_garch_model)$ll3)
    # aic
    AIC1 <- -2 * logLik(gnar_model)$ll1 + 2 * attr(logLik(gnar_model),'df')
    AIC2 <- -2 * logLik(gnar_garch_model)$ll3 + 2 * attr(logLik(gnar_garch_model),'df')
    # parameters
    parameters1 <- gnar_model$mod$coef
    parameters2 <- gnar_garch_model$mod$coef
    # matrix
    matrix[i,] <- c(logLiks1, logLiks2, AIC1, AIC2, parameters1, parameters2)
}

# mean of columns
results_mean <- round(colMeans(matrix, na.rm=TRUE),4)
list(logLiks1=results_mean[1:3], logLiks2=results_mean[4:6], AIC1=results_mean[7], AIC2=results_mean[8], 
    parameters1=results_mean[9:(8+length(global_alpha)+sum(lengths(beta)))], 
    parameters2=results_mean[(9+length(global_alpha)+sum(lengths(beta))):(8+2*length(global_alpha)+2*sum(lengths(beta)))])

# std of columns
results_std <- round(apply(matrix, 2, sd, na.rm=TRUE),4)
list(logLiks1=results_std[1:3], logLiks2=results_std[4:6], AIC1=results_std[7], AIC2=results_std[8], 
    parameters1=results_std[9:(8+length(global_alpha)+sum(lengths(beta)))], 
    parameters2=results_std[(9+length(global_alpha)+sum(lengths(beta))):(8+2*length(global_alpha)+2*sum(lengths(beta)))])