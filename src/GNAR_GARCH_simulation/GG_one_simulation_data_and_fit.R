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

# Plot the network
par(mar = c(1, 1, 1, 1))
plot(fiveNet, vertex.label.cex=1.5)
dev.copy2pdf(file=here('figures','GNAR_GARCH_simulation','GG_network_plot.pdf'), width=8, height=8)
par(mar = c(5, 4, 4, 2) + 0.1)

# Parameters
global_alpha <- c(0.3, 0.2)
beta <- list(c(0.1, 0.05),c(0.1))
mu <- rep(0.1,5)
a <- list(c(0.5),c(0.4),c(0.3),c(0.3),c(0.4))
b <- list(c(0.4),c(0.3),c(0.2),c(0.3),c(0.4))

# Simulation Data
set.seed(11111)
data <- simulate_gnar_garch(n=n, global_alpha=global_alpha, beta=beta, mu=mu, a=a, b=b)

# Fit a GNAR model
gnar_model <- GNARfit(vts=data, net=fiveNet, alphaOrder=length(global_alpha),betaOrder=lengths(beta))
# Fit a GNAR + GARCH model
gnar_garch_model <- GNARfit(vts=data, net=fiveNet, alphaOrder=length(global_alpha),betaOrder=lengths(beta),archOrder=length(a[1]),garchOrder=length(b[1]))

# Fitted values and std of residuals for GNAR model
fitted1 <- ts(fitted(gnar_model)[,1], end=n, frequency=1)
# Fitted values and conditional std for GNAR + GARCH model
fitted2 <- ts(fitted(gnar_garch_model)[,1], end=n, frequency=1)
spec <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(length(a[1]), length(b[1]))),
                    mean.model = list(armaOrder = c(length(global_alpha), 0),include.mean=FALSE),
                    distribution = 'norm')
fit <- ugarchfit(spec = spec, data = residuals(gnar_model)[,1])
conditional_std <- ts(sigma(fit), end=n, frequency=1)

# Save as cvs
write.csv(data, file=here('data','derived','GNAR_GARCH_simulation','GG_one_simulation_data.csv'))
write.csv(fitted(gnar_model), file=here('data','derived','GNAR_GARCH_simulation','GG_one_simulation_fittedvalues_gnar.csv'))
write.csv(fitted(gnar_garch_model), file=here('data','derived','GNAR_GARCH_simulation','GG_one_simulation_fittedvalues_gnar_garch.csv'))
write.csv(residuals(gnar_model), file=here('data','derived','GNAR_GARCH_simulation','GG_one_simulation_residuals_gnar.csv'))
write.csv(conditional_std, file=here('data','derived','GNAR_GARCH_simulation','GG_one_simulation_conditional_std_gnar_garch_node1.csv'))