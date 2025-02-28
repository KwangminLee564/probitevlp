#mvprobit.R


#'
#' @export
#'
#' @examples
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' inputdata <- generate_data(param,2000)
#' mcmc.out <- mvprobit(inputdata)
#'
mvprobit <- function(inputdata,postnum=5000){

  library(magrittr)
  r <- dim(inputdata$Y)[2]
  p <- dim(inputdata$X)[2]+1

  Data1 <- list(p=r , X=cbind(1,inputdata$X) %x% diag(1,r),y=as.vector(t(inputdata$Y)))
  Mcmc1 = list(R=postnum, keep=1)

  out = bayesm::rmvpGibbs(Data=Data1, Mcmc=Mcmc1)


  Sigsample <- 1:dim(out$sigmadraw)[1] %>%
    purrr::map(~matrix(out$sigmadraw[.x,],r,r)) %>%
    abind::abind(along=3)

  Betasample <- 1:dim(out$sigmadraw)[1] %>%
    purrr::map(~matrix(out$betadraw[.x,],r,p)) %>%
    abind::abind(along=3)

  list(Betasample= Betasample[,,-(1:as.integer(postnum/2))],
       Sigsample=Sigsample[,,-(1:as.integer(postnum/2))])
}





