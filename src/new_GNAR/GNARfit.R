GNARfit <- function(vts=GNAR::fiveVTS, net=GNAR::fiveNet, alphaOrder=2, betaOrder=c(1,1), archOrder=0, garchOrder=0, fact.var=NULL,
                    globalalpha=TRUE, tvnets=NULL, netsstart=NULL, ErrorIfNoNei=TRUE){
  #input checks
  stopifnot(is.GNARnet(net))
  stopifnot(ncol(vts) == length(net$edges))
  stopifnot(alphaOrder > 0)
  stopifnot(floor(alphaOrder) == alphaOrder)
  stopifnot(length(betaOrder) == alphaOrder)
  stopifnot(floor(betaOrder) == betaOrder)
  if(!is.null(fact.var)){
    stopifnot(length(fact.var) == length(net$edges))
  }
  stopifnot(is.matrix(vts))
  stopifnot(is.logical(globalalpha))
  if(!is.null(tvnets)){
    cat("Time-varying networks not yet supported")
  }
  stopifnot(is.null(tvnets))
  useNofNei <- 1
  # end of input checks
  frbic <- list(nnodes=length(net$edges),alphas.in=alphaOrder,betas.in=betaOrder,arch.in=archOrder,garch.in=garchOrder,fact.var=fact.var,
                globalalpha=globalalpha,xtsp=tsp(vts),time.in=nrow(vts),net.in=net,final.in=vts[(nrow(vts)-alphaOrder+1):nrow(vts),])
  # design matrix X
  dmat <- GNARdesign(vts=vts, net=net, alphaOrder=alphaOrder, betaOrder=betaOrder, fact.var=fact.var,
                     globalalpha=globalalpha, tvnets=tvnets, netsstart=netsstart)
  if(ErrorIfNoNei){
    if(any(apply(dmat==0, 2, all))){
      parname <- strsplit(names(which(apply(dmat==0, 2, all)))[1], split=NULL)[[1]]
      betastage <- parname[(which(parname==".")+1) :(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", betastage)
    }
  }
  predt <- nrow(vts)-alphaOrder # T-p
  # observation vector y and weights vector
  yvec <- NULL
  mu <- NULL
  arch.coefs <- NULL
  garch.coefs <- NULL
  sq_sigma <- NULL
  weights <- NULL
  for(ii in 1:length(net$edges)){
    yts <- vts[((alphaOrder+1):(predt+alphaOrder)),ii]
    yvec <- c(yvec, yts)
    if (archOrder>0 || garchOrder>0){
      spec <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(archOrder, garchOrder)),
                        mean.model = list(armaOrder = c(alphaOrder, 0),include.mean=FALSE),
                        distribution = 'norm')
      fit <- ugarchfit(spec = spec, data = yts)
      # get the coefficients
      mu <- c(mu, coef(fit)[alphaOrder+1])
      if (archOrder != 0) {
        arch.coefs <- c(arch.coefs, coef(fit)[(alphaOrder+2):(alphaOrder+1+archOrder)])
      }
      if (garchOrder != 0) {
        garch.coefs <- c(garch.coefs, coef(fit)[(alphaOrder+2+archOrder):(alphaOrder+1+archOrder+garchOrder)])
      }
      # get the conditional variances and residuals
      sq_sigma <- c(sq_sigma, as.vector(sigma(fit)^2))
      weights <- c(weights, 1/as.vector(sigma(fit)^2))
    }
  }
  # remove NA
  if(sum(is.na(yvec))>0){
    yvec2 <- yvec[!is.na(yvec)]
    dmat2 <- dmat[!is.na(yvec),]
    modNoIntercept <- lm(yvec2~dmat2+0)
  }else{
    modNoIntercept <- lm(yvec~dmat+0, weights=weights) # fit model
  }
  out <- list(mod=modNoIntercept, y=yvec, dd=dmat, weights=weights, 
              sq_sigma=sq_sigma, mu=mu, arch.coefs=arch.coefs, 
              garch.coefs=garch.coefs, frbic=frbic)
  class(out) <- "GNARfit"
  return(out)
}