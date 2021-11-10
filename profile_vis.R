library(profvis)

profvis({
  r <- 3
  p <- 4
  u <- u.tru <- 2
  param <- generate_par(r, p, u)
  inputdata <- generate_data(param,100)
  mcmc.env <- envprobit_BIC(inputdata,u,mcmc.num=10)

})
