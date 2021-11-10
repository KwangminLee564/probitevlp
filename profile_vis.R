library(profvis)

profvis({
  r <- 20
  p <- 10
  u <- u.tru <- 2
  n <- 200
  param <- generate_par(r, p, u)
  inputdata <- generate_data(param,n)
  mcmc.env <- envprobit_BIC(inputdata,u,mcmc.num=10000)

})
