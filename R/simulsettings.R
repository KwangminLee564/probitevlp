

#'
#'
#' @export
#'
#' @examples
#'
#' r <- 20
#' p <- 10
#' u <- u.tru <- 2
#' all_pars <- generate_par(r, p, u)
#'
#'
generate_par <- function(r, p, u) {
  mu.tru <- runif(r, -1, 1)
  eta.tru <- matrix(runif(u * p, min = 0, max = 1), nrow = u, ncol = p)
  A <- A.tru <- matrix(runif(u * (r - u), min = -1, max = 1), nrow = r-u, ncol = u)
  gamma_gamma0 <- find_gammas_from_A(A)
  gamma <- gamma_gamma0$gamma
  gamma0 <- gamma_gamma0$gamma0
  CA <- gamma_gamma0$CA

  # CA <- rbind(diag(1, u), A)
  # DA <- rbind(-t(A), diag(1, r-u))
  #
  # calculate beta
  beta.tru <- CA %*% eta.tru
  #
  # CAtCA <- t(CA) %*% CA
  # gamma <- CA %*% sqrtmatinv(CAtCA)
  #
  # DAtDA <- t(DA) %*% DA
  # gamma0 <- DA %*% sqrtmatinv(DAtDA)

  # omega.tru <- runif(0.5, 0.8)
  # Omega.tru <- rinvwish(dim = u, Phi = diag(0.8, u), nu = 1)
  # Omega0.tru <- rinvwish(dim = r-u, Phi = diag(250, r-u), nu = 1)
  # Omega.tru <- diag(sort(runif(u, 0, 1), decreasing = TRUE),
  #                   ncol = u, nrow = u)
  # Omega0.tru <- diag(sort(runif(r-u, 2, 10), decreasing = TRUE),
  #                    ncol = r-u, nrow = r-u)

  Omega.tru <- diag(sort(runif(u, 0, 1), decreasing = TRUE),
                    ncol = u, nrow = u)
  Omega0.tru <- diag(sort(runif(r-u, 5, 10), decreasing = TRUE),
                     ncol = r-u, nrow = r-u)


  Sigma1 <- gamma %*% Omega.tru %*% t(gamma)
  Sigma2 <- gamma0 %*% Omega0.tru %*% t(gamma0)

  Sigma.tru <- Sigma1 + Sigma2

  list(mu.tru = mu.tru,
       eta.tru=eta.tru,
       beta.tru = beta.tru,
       Omega.tru = Omega.tru,
       Omega0.tru = Omega0.tru,
       A.tru = A.tru,
       gamma.tru = gamma,
       gamma0.tru = gamma0,
       Sigma.tru = Sigma.tru,
       mux = 0,
       sigmax = 1)
}


#'
#' @export
#'
#' @examples
#'
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' param <- generate_par(r, p, u)
#' simuldata <- generate_data(param,50)
#'
#'
generate_data <- function(param,n,  ...)
{
  mux <- param$mux
  sigmax <- param$sigmax
  mu.tru <- param$mu.tru
  Sigma.tru <- param$Sigma.tru
  beta.tru <- param$beta.tru
  p <- dim(beta.tru)[2]
  r <- dim(beta.tru)[1]



  X <- matrix(rnorm(n*p, mean = mux, sd = sigmax),
              nrow = n, ncol = p)
  eps <- matrix(rnorm(n*r), n, r) %*% sqrtmat(Sigma.tru)
  Z <- tcrossprod(rep(1, n), mu.tru) + X %*% t(beta.tru) + eps
  Y <- apply(Z>0,2,as.integer)

  list(X = X, Y = Y,Z=Z)
}

