logLik.GNARfit <- function(object,...){
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  t <- object$frbic$time.in
  n <- object$frbic$nnodes
  stopifnot(floor(t)==t)
  stopifnot(floor(n)==n)

  # eps <- residuals(object)
  # eps[is.na(eps)] <- 0
  # Sigma <- crossprod(eps)/t # t(eps) %*% eps / t
  # ll <- -(t*n) * log(2*pi)/2 - (t/2)*log(det(Sigma)) - sum(diag(eps %*% solve(Sigma) %*% t(eps)))/2

  # error_vector <- object$mod$residuals
  # if (is.null(object$mod$weights)) {
  #   var <- sum(error_vector^2)/(t*n)
  #   ll <- -(t*n)*log(2*pi)/2 -(t*n)*log(var)/2 - (t*n)/2
  # } else {
  #   sq_sigma <- object$sq_sigma
  #   ll <- -(t*n)*log(2*pi)/2 - sum(log(sq_sigma))/2 - sum(error_vector^2/sq_sigma)/2
  # }

  error_vector <- object$mod$residuals
  eps <- residuals(object)
  eps[is.na(eps)] <- 0
  Sigma <- crossprod(eps)/t # t(eps) %*% eps / t
  var <- sum(error_vector^2)/(t*n)
  if (is.null(object$mod$weights)) {
    ll1 <- -(t*n)*log(2*pi)/2 -(t*n)*log(var)/2 - (t*n)/2
    ll2 <- -(t*n) * log(2*pi)/2 - (t/2)*log(det(Sigma)) - (t*n)/2
    ll3 <- NA
  } else {
    ll1 <- -(t*n)*log(2*pi)/2 -(t*n)*log(var)/2 - (t*n)/2
    ll2 <- -(t*n) * log(2*pi)/2 - (t/2)*log(det(Sigma)) - (t*n)/2
    sq_sigma <- object$sq_sigma
    ll3 <- -(t*n)*log(2*pi)/2 - sum(log(sq_sigma))/2 - sum(error_vector^2/sq_sigma)/2
  }

  lls <- list(ll1=ll1,ll2=ll2,ll3=ll3)
  class(lls) <- "logLiks"
  attr(lls, "df") <- length(coef(object)) + length(object$mu) + length(object$arch.coefs) + length(object$garch.coefs)
  return(lls)
}