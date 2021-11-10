
#'
#' @export
#' @examples
#'
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' all_pars <- generate_par(r, p, u)
#' simuldata <- generate_data(all_pars,50)
#' A0 <- matrix(0,r-u,u)
#' K <- diag(10^4, r-u)
#' L <- diag(10^4, u)
#' M <- diag(0.001,p)
#' K.half.inv <- sqrtmatinv(K)
#' L.half.inv <- sqrtmatinv(L)
#' M.half <- sqrtmat(M)
#'
#' Wmat <- simuldata$Z
#' Xmat <- simuldata$X
#' A <- all_pars$A.tru
#' teta <- all_pars$eta.tru
#' tmu <- all_pars$mu.tru
#' Omega.inv <- solve(all_pars$Omega.tru)
#' Omega0.inv <- solve(all_pars$Omega0.tru)
#' lpd_A_pred(A,Wmat,Xmat,teta,tmu,Omega0.inv,Omega.inv,K.half.inv, L.half.inv, A0,M.half)
#'
lpd_A_pred <- function(A, Wmat, Xmat,
                       teta,tmu,
                       Omega0.inv,
                       Omega.inv,
                       K.half.inv, L.half.inv, A0,
                       M.half) {

  p <- dim(Xmat)[2]
  gamma_gamma0 <- find_gammas_from_A(A)
  gamma <- gamma_gamma0$gamma
  gamma0 <- gamma_gamma0$gamma0

  term1 <- -p*as.double(determinant(gamma_gamma0$CAtCA)$modulus)
  term2 <- sum((K.half.inv %*% (A-A0) %*% L.half.inv)^2)

  Sigma.inv <- gamma %*% tcrossprod(Omega.inv, gamma) +
    gamma0 %*% tcrossprod(Omega0.inv, gamma0)
  resi <- Wmat - t(as.vector(tmu) + gamma_gamma0$CA %*% teta %*% t(Xmat))
  term3 <-  sum((resi %*% Sigma.inv) * resi)

  term4 <- sum(Omega.inv * crossprod(M.half %*% t(teta) %*% sqrtmatinv(gamma_gamma0$CAtCA)))
  -0.5 * (term1 + term2 + term3 + term4)
}


#'
#' @examples
#'
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' all_pars <- generate_par(r, p, u)
#' simuldata <- generate_data(all_pars,50)
#' A0 <- matrix(0,r-u,u)
#' K <- diag(10^4, r-u)
#' L <- diag(10^4, u)
#' M <- diag(0.001,p)
#' K.half.inv <- sqrtmatinv(K)
#' L.half.inv <- sqrtmatinv(L)
#' M.half <- sqrtmat(M)
#'
#' Wmat <- simuldata$Z
#' Xmat <- simuldata$X
#' A <- all_pars$A.tru
#' teta <- all_pars$eta.tru
#' tmu <- all_pars$mu.tru
#' Omega.inv <- solve(all_pars$Omega.tru)
#' Omega0.inv <- solve(all_pars$Omega0.tru)
#' lpd_val <- lpd_A_pred(A,Wmat,Xmat,teta,tmu,Omega0.inv,Omega.inv,K.half.inv, L.half.inv, A0,M.half)
#'
#' rw_var <- rep(0.4,dim(A)[2])
#' rw_var <- c(0.6,0.3)
#' rwmh_colwise(lpd_val, rw_var,A,Wmat, Xmat,teta,tmu,Omega0.inv,
#' Omega.inv,K.half.inv, L.half.inv, A0,M.half)
#'
#'
rwmh_colwise <- function(lpd_val, rw_var,A,...){

  #if(is.null(lpd_val)){
  #  lpd_val <- lpd_fun(A=A,...)
  #}
  u <- dim(A)[2]
  r_u <- dim(A)[1]
  alphas <- numeric(u)

  for(j in sample(1:u)) {
    # if (autotune_size == "single") tau_curr <- tau
    A_j_star <- A[, j] + rnorm(r_u, 0, sqrt(rw_var[j]))
    A_star <- A
    A_star[, j] <- A_j_star
    lpd_A_star <- lpd_A_pred(A = A_star, ...)
    alphas[j] <- min(exp(lpd_A_star-lpd_val),1)
    if (runif(1) < alphas[j]) {
      # accept
      A <- A_star
      lpd_val <- lpd_A_star
    }
  }

  list(A=A,lpd_val=lpd_val,alphas=alphas)
}

