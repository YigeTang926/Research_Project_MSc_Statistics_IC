library(here)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
i_am("src/data_analyses/TBATS_ average_daily.R")
library(igraph)
library(GNAR)
library(rugarch)
library(ggplot2)
library(viridis)
library(plotly)

source(here('src','new_GNAR','GNARfit.R'))
source(here('src','new_GNAR','logLik.GNARfit.R'))

# Load data
network_nodes <- read.csv(file=here('data','raw','network_nodes.csv'))
network_edges <- read.csv(file=here('data','raw','network_edges.csv'))
wind_remainder <- read.csv(file=here('data','derived','czech','seasonal_decomposition','wind_remainder.csv') ,row.names=1)
solar_remainder <- read.csv(file=here('data','derived','czech','seasonal_decomposition','solar_remainder.csv') ,row.names=1)

# Subset the network nodes and edges in Czech
network_nodes_czech <- subset(network_nodes, country == "CZE")
network_edges_czech <- network_edges[network_edges[,"fromNode"] %in% network_nodes_czech[,'ID'] & network_edges[,"toNode"] %in% network_nodes_czech[,'ID'],][,c("fromNode","toNode","length")]

# Create a graph object
graph <- graph_from_data_frame(network_edges_czech, directed=TRUE, vertices=network_nodes_czech[,'ID'])
# Add the weights, which is 1/length
E(graph)$weight <- 1/E(graph)$length
# Plot the graph
layout <- layout_nicely(graph)
plot(graph, vertex.label=V(graph)$name, vertex.size=11, vertex.label.cex=0.8, edge.arrow.size=0.6, 
    edge.label.cex=0.8, edge.label=round(E(graph)$length,0), edge.width=E(graph)$weight*120, layout=layout)
dev.copy2pdf(file=here('figures','czech','igraph_czech.pdf'), width=8, height=8)

# Convert the igraph object into a GNAR network object
network <- igraphtoGNAR(graph)

# Convert the remainder into a ts
wind_remainder_ts <- ts(wind_remainder, end=c(2014,365), frequency=365)
solar_remainder_ts <- ts(solar_remainder, end=c(2014,365), frequency=365)

# GNAR(2,{2,1}) + GARCH(1,1) Model
alpha_order <- 2
v <- 1
u <- 1
wind_model1 <- GNARfit(vts=wind_remainder_ts, net=network, alphaOrder=alpha_order, betaOrder=c(2,1), globalalpha=TRUE)
wind_model2 <- GNARfit(vts=wind_remainder_ts, net=network, alphaOrder=alpha_order, betaOrder=c(2,1), archOrder=v, garchOrder=u, globalalpha=TRUE)
solar_model1 <- GNARfit(vts=solar_remainder_ts, net=network, alphaOrder=alpha_order, betaOrder=c(2,1), globalalpha=TRUE)
solar_model2 <- GNARfit(vts=solar_remainder_ts, net=network, alphaOrder=alpha_order, betaOrder=c(2,1), archOrder=v, garchOrder=u, globalalpha=TRUE)

# Get the log-likelihood and AIC
logLik_AIC <- matrix(0, nrow=4, ncol=4)
logLik_AIC[1,] <- c(logLik(wind_model1)$ll1, logLik(wind_model1)$ll2, logLik(wind_model1)$ll3, -2 * logLik(wind_model1)$ll1 + 2 * attr(logLik(wind_model1),'df'))
logLik_AIC[2,] <- c(logLik(wind_model2)$ll1, logLik(wind_model2)$ll2, logLik(wind_model2)$ll3, -2 * logLik(wind_model2)$ll3 + 2 * attr(logLik(wind_model2),'df'))
logLik_AIC[3,] <- c(logLik(solar_model1)$ll1, logLik(solar_model1)$ll2, logLik(solar_model1)$ll3, -2 * logLik(solar_model1)$ll1 + 2 * attr(logLik(solar_model1),'df'))
logLik_AIC[4,] <- c(logLik(solar_model2)$ll1, logLik(solar_model2)$ll2, logLik(solar_model2)$ll3, -2 * logLik(solar_model2)$ll3 + 2 * attr(logLik(solar_model2),'df'))
colnames(logLik_AIC) <- c('logLik1','logLik2','logLik3','AIC')
write.csv(logLik_AIC, file=here('data','derived','czech','gnar_garch','gnar_garch_logLik_AIC.csv'))

# Get the individual log-likelihood
n_nodes <- length(V(graph)$name)
T <- length(wind_remainder_ts[,1])
individual_logLik <- matrix(0, nrow=n_nodes, ncol=4)
for (i in 1:n_nodes) {
    logLik_gnar <- -T*log(2*pi)/2 -T*log(sum(wind_model1$mod$residuals^2)/(n_nodes*T))/2 - T/2
    sq_sigma <- wind_model2$sq_sigma[((i-1)*(T-alpha_order)+1):(i*(T-alpha_order))]
    logLik_gnar_garch <- -T*log(2*pi)/2 - sum(log(sq_sigma))/2 - sum(residuals(wind_model2)[,i]^2/sq_sigma)/2
    individual_logLik[i,1] <- logLik_gnar
    individual_logLik[i,2] <- logLik_gnar_garch
    logLik_gnar <- -T*log(2*pi)/2 -T*log(sum(solar_model1$mod$residuals^2)/(n_nodes*T))/2 - T/2
    sq_sigma <- solar_model2$sq_sigma[((i-1)*(T-alpha_order)+1):(i*(T-alpha_order))]
    logLik_gnar_garch <- -T*log(2*pi)/2 - sum(log(sq_sigma))/2 - sum(residuals(solar_model2)[,i]^2/sq_sigma)/2
    individual_logLik[i,3] <- logLik_gnar
    individual_logLik[i,4] <- logLik_gnar_garch
}
write.csv(individual_logLik, file=here('data','derived','czech','gnar_garch','individual_logLik.csv'))

# Node 1061
GG_results_1061 <- matrix(0, nrow=T, ncol=6)
GG_results_1061[,1] <- wind_remainder_ts[,22]
GG_results_1061[-(1:alpha_order),2] <- fitted(wind_model2)[,22]
GG_results_1061[-(1:alpha_order),3] <- wind_model2$sq_sigma[((22-1)*(T-alpha_order)+1):(22*(T-alpha_order))]
GG_results_1061[,4] <- solar_remainder_ts[,22]
GG_results_1061[-(1:alpha_order),5] <- fitted(solar_model2)[,22]
GG_results_1061[-(1:alpha_order),6] <- solar_model2$sq_sigma[((22-1)*(T-alpha_order)+1):(22*(T-alpha_order))]
colnames(GG_results_1061) <- c('observations','fitted values','sq_sigma','observations','fitted values','sq_sigma')
write.csv(GG_results_1061, file=here('data','derived','czech','gnar_garch','GG_results_1061.csv'))

# Select alpha and beta using AIC 
l <- 10
m <- 7
wind_AIC_1 <- matrix(0, nrow = l, ncol = m)
wind_AIC_2 <- matrix(0, nrow = l, ncol = m)
for(a in 1:l) {
    for(b in 0:(m-1)) {
        model <- GNARfit(vts=wind_remainder_ts,  net=network, alphaOrder=a, betaOrder=rep(b,a), archOrder=1, garchOrder=1, globalalpha=TRUE)
        wind_AIC_1[a, b + 1] <- -2 * logLik(model)$ll2 + 2 * attr(logLik(model), 'df')
        wind_AIC_2[a, b + 1] <- -2 * logLik(model)$ll3 + 2 * attr(logLik(model), 'df')
    }    
}
solar_AIC_1 <- matrix(0, nrow = l, ncol = m)
solar_AIC_2 <- matrix(0, nrow = l, ncol = m)
for(a in 1:l) {
    for(b in 0:(m-1)) {
        model <- GNARfit(vts=solar_remainder_ts,  net=network, alphaOrder=a, betaOrder=rep(b,a), archOrder=1, garchOrder=1, globalalpha=TRUE)
        solar_AIC_1[a, b + 1] <- -2 * logLik(model)$ll2 + 2 * attr(logLik(model), 'df')
        solar_AIC_2[a, b + 1] <- -2 * logLik(model)$ll3 + 2 * attr(logLik(model), 'df')
    }    
}
# Save the AIC
write.csv(wind_AIC_1, file=here('data','derived','czech','gnar_garch_AIC','wind_AIC_1.csv'))
write.csv(wind_AIC_2, file=here('data','derived','czech','gnar_garch_AIC','wind_AIC_2.csv'))
write.csv(solar_AIC_1, file=here('data','derived','czech','gnar_garch_AIC','solar_AIC_1.csv'))
write.csv(solar_AIC_2, file=here('data','derived','czech','gnar_garch_AIC','solar_AIC_2.csv'))

# Select beta when alpha=2 using AIC
m <- 7
wind_AIC_beta_1 <- matrix(0, ncol = m, nrow = m)
wind_AIC_beta_2 <- matrix(0, ncol = m, nrow = m)
for(b1 in 0:(m-1)) {
    for(b2 in 0:(m-1)) {
        model <- GNARfit(vts=wind_remainder_ts,  net=network, alphaOrder=2, betaOrder=c(b1, b2), archOrder=1, garchOrder=1, globalalpha=TRUE)
        wind_AIC_beta_1[b1 + 1, b2 + 1] <- -2 * logLik(model)$ll2 + 2 * attr(logLik(model), 'df')
        wind_AIC_beta_2[b1 + 1, b2 + 1] <- -2 * logLik(model)$ll3 + 2 * attr(logLik(model), 'df')
    }    
}
solar_AIC_beta_1 <- matrix(0, ncol = m, nrow = m)
solar_AIC_beta_2 <- matrix(0, ncol = m, nrow = m)
for(b1 in 0:(m-1)) {
    for(b2 in 0:(m-1)) {
        model <- GNARfit(vts=solar_remainder_ts,  net=network, alphaOrder=2, betaOrder=c(b1, b2), archOrder=1, garchOrder=1, globalalpha=TRUE)
        solar_AIC_beta_1[b1 + 1, b2 + 1] <- -2 * logLik(model)$ll2 + 2 * attr(logLik(model), 'df')
        solar_AIC_beta_2[b1 + 1, b2 + 1] <- -2 * logLik(model)$ll3 + 2 * attr(logLik(model), 'df')
    }    
}

# Save the AIC
write.csv(wind_AIC_beta_1, file=here('data','derived','czech','gnar_garch','wind_AIC_beta_1.csv'))
write.csv(wind_AIC_beta_2, file=here('data','derived','czech','gnar_garch','wind_AIC_beta_2.csv'))
write.csv(solar_AIC_beta_1, file=here('data','derived','czech','gnar_garch','solar_AIC_beta_1.csv'))
write.csv(solar_AIC_beta_2, file=here('data','derived','czech','gnar_garch','solar_AIC_beta_2.csv'))