#' @export
fitENIGMA <- function(Y,formula=NULL,X=NULL,L=2,alpha=NULL,data=NULL,...){
  if(is.null(X)){
    X <- model.matrix(formula,data = data)
  }
  if(is.null(alpha)){
    alpha=rep(1,L)
  }
  dat4stan <- list(N=nrow(Y),K=ncol(Y),D=ncol(X),
                   L=L,
                   y=Y,
                   X=X,
                   alpha=alpha)
  rstan::optimizing(stanmodels$softmaxreg,dat4stan,as_vector=FALSE,hessian=TRUE,...)
}

#' @export
calcSE <- function(fit){
  se <-sqrt(-diag(solve(fit$hessian)))
  df <-data.frame(param=names(se),se)
  rownames(df) <- NULL
  return(df)
}

#' @export
calcEvidence <- function(fit){
  d0 <-length(fit$par$pi) +
    length(fit$par$beta) +
    length(fit$par$gamma) +
    length(fit$par$sigma) +
    length(fit$par$tau)
  ldet0 <- determinant(fit$hessian)
  fit$value+d0*log(2*pi)-ldet0$sign*ldet0$modulus[1]/2
}
