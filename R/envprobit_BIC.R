#'
#' @export
#'
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,100)
#' mcmc.env <- envprobit_BIC(inputdata,u,mcmc.num=20)
#'
envprobit_BIC <- function(inputdata,u,...){
  mcmc.out <- envprobit(inputdata,u,...)
  1:dim(mcmc.out$Betasample)[3] %>%
    purrr::map(function(ind){
      loglik_probit(mcmc.out$Betasample[,,ind],
                    mcmc.out$Sigsample[,,ind],
                    inputdata)}) %>% unlist

  #xx = maximized likelihood
  # -2* xx + (r*(r + 1)/2 + r + p*u-1)*log(n)


  #if(ncore>1){
  #  library(furrr)
  #  plan(multisession, workers = ncore)
  #  furrr::future_map(uvec,~envprobit(inputdata,.x,...))
  #} else{
  #  purrr::map(uvec,~envprobit(inputdata,.x,...))
  #}


}



#'
#' @export
#'
#'
#'
calc_BIC <- function(mubeta,Sigma,inputdata,u){
  p <- dim(inputdata$X)[2]
  r <- dim(inputdata$Y)[2]
  n <- dim(inputdata$X)[1]
  -2*loglik_probit(mubeta,Sigma,inputdata) + (r*(r + 1)/2 + r + p*u)*log(n)
}


#'
#' @export
#'
#' @examples
#'
#' r <- 3
#' p <- 5
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,50)
#' loglik_probit(cbind(param$mu.tru,param$beta.tru),param$Sigma.tru,inputdata)
#'
loglik_probit <- function(mubeta,Sigma,inputdata){
  p <- dim(inputdata$X)[2]
  r <- dim(inputdata$Y)[2]

  colnames(inputdata$X) <- paste0("x",1:p)
  xMat <- cbind(const=1,inputdata$X)

  yMat <- inputdata$Y
  colnames(yMat) <- paste0("y",1:r)
  formula_str <- paste0("cbind(",paste(colnames(yMat),collapse = ","),")",
                        "~",paste(paste0("x",1:p),collapse = "+"))

  Sigma <- LaplacesDemon::as.positive.definite(Sigma)
  sum(mvProbit::mvProbitLogLik(as.formula(formula_str),data=as.data.frame(cbind(xMat,yMat)),
                               coef=c(t(scaletoCor(mubeta,Sigma))),sigma=cov2cor(Sigma)))
}

