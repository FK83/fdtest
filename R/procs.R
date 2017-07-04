#' @import murphydiagram
#' @export fdtest
#' @export sim_dgp
#' @export sim_dgp_est
#' @export sim_dgp_quantiles
#' @importFrom stats arima.sim rnorm runif qnorm lm sd
#' @useDynLib fdtest

# Sign helper function
signhelper <- function(mat, signs){
  mat[signs == -1, ] <- -mat[signs == -1, ]
  mat
}

#' Randomization test
#'
#'@description Randomization test for forecast dominance
#'@param fc matrix with forecasts (first two cols) and realization (third col)
#'@param np nr of randomization draws
#'@param functional either "expectile" or "quantile"
#'@param alpha level of the expectile or quantile
#'@param test_functional test statistic to be used, either "sqpos" or "abspos"
#'@details \emph{Note:} The test's H0 is that forecast 1 weakly dominates forecast 2. Hence, low p-values indicate that forecast 1 is \emph{not} dominant. This can mean that either i) forecast 2 is dominant or that ii) none of the forecasts dominates the other.
fdtest <- function(fc, functional = "expectile", alpha = 0.5,
                   np = 1000, test_functional = "sqpos"){

  # Starting time
  t0 <- Sys.time()

  # Input check
  if (!is.matrix(fc)){
    fc <- as.matrix(fc)
  }

  # Functional of interest
  S <- function(x, y, theta){
    sapply(theta,
           function(z) extremal_score(x, y, theta = z,
                                      functional = functional,
                                      alpha = alpha))
  }

  # Grid for theta
  theta_grid <- seq(from = min(fc) - 0.1 * sd(fc), to = max(fc) +
                0.1 * sd(fc), length.out = 100)

  # Sampe size
  T <- nrow(fc)

  # Compute score differences
  S_diff <- S(x = fc[, 1], y = fc[, 3], theta = theta_grid) -
    S(x = fc[, 2], y = fc[, 3], theta = theta_grid)

  # Function for test stat
  if (test_functional == "abspos"){
    teststat <- function(z){
      sum(z[z >= 0])
    }
    auxin <- 2
  } else if (test_functional == "sqpos"){
    teststat <- function(z){
      sum(z[z >= 0]^2)
    }
    auxin <- 1
  }

  # Compute test statistic (functional of theta)
  teststat_sample <- teststat(colMeans(S_diff))

  # Loop
  Bpy <- permhelper(S_diff, auxin, np)

  # p-value
  pval <- mean(Bpy > teststat_sample)

  # timing
  t1 <- Sys.time()

  # Return list
  list(teststat = teststat_sample, pval = pval, time = (t1 - t0))
}

#' Functions to simulate example data
#'
#' @param n nr of observations
#' @param g persistence parameter
#' @param sd_2 variance parameter, determines relative performance
#' @param h forecast horizon
#' @return matrix of dimension (n x 3), containing forecasts (columns 1, 2) and realizations (column 3) from example process
#' @rdname sim_dgp
sim_dgp <- function(n, g = 0.8, h = 1, sd_2 = 1){
  x <- matrix(0, n, 2)
  if (g != 0){
    x[, 1] <- arima.sim(model = list(ar = g), n = n)
    x[, 2] <- arima.sim(model = list(ar = g), n = n, sd = sd_2)
  } else {
    x[, 1] <- rnorm(n)
    x[, 2] <- rnorm(n, sd = sd_2)
  }
  f <- matrix(0, n - h, 2)
  f[, 1] <- x[(h+1):n, 1] + x[1:(n-h), 2]*(g^h)
  f[, 2] <- x[(h+1):n, 2] + x[1:(n-h), 1]*(g^h)
  y <- apply(x, 1, sum)
  return(cbind(f, y[-(1:h)]))
}

#' @param presample, number of presample observations
#' @rdname sim_dgp
sim_dgp_est <- function(n, g = 0.8, sd_2 = 1, presample = 50){
  n2 <- n + presample
  x <- matrix(0, n2, 2)
  x[, 1] <- arima.sim(model = list(ar = g), n = n2)
  x[, 2] <- arima.sim(model = list(ar = g), n = n2, sd = sd_2)
  f <- matrix(0, n, 2)
  for (ii in presample:(n2-1)){
    aux <- cbind(c(x[2:ii, ]), c(x[1:(ii-1), ]))
    ghat <- sum(aux[,1]*aux[,2])/sum(aux[,2]^2)
    f[ii - presample + 1, 1] <- x[ii + 1, 1] + x[ii, 2]*ghat
    f[ii - presample + 1, 2] <- x[ii + 1, 2] + x[ii, 1]*ghat
  }
  y <- apply(x, 1, sum)
  return(cbind(f, y[-(1:presample)]))
}

sim_ar1garch11 <- function(n, bi = 500){

  n2 <- n + bi
  e <- rnorm(n2)
  s2 <- rep(1, n2)
  y <- rep(0, n2)
  m <- rep(NA, n2)
  for (tt in 2:n2){
    m[tt] <- 0.03 + 0.05*y[tt-1]
    s2[tt] <- 0.05 + 0.9*s2[tt-1] + 0.05*(e[tt-1]^2)*s2[tt-1]
    y[tt] <- m[tt] + sqrt(s2[tt])*e[tt]
  }
  data.frame(m = m[-(1:bi)], s2 = s2[-(1:bi)], y = y[-(1:bi)])

}

#' @param a quantile level
#' @param tau1 precision parameter for forecaster 1
#' @param tau2 precision parameter for forecaster 2
#' @aliases sim_dgp_quantiles, sim_dgp_est
#' @rdname sim_dgp
sim_dgp_quantiles <- function(n, a = 0.05, tau1, tau2){

  dat <- sim_ar1garch11(n)
  mat <- matrix(0, n, 3)
  mat[, 3] <- dat$y
  aux1 <- log(a) - log(1-a) + rnorm(n, sd = tau1)
  mat[, 1] <- qnorm(1/(1+exp(-aux1)), mean = dat$m, sd = sqrt(dat$s2))
  aux2 <- log(a) - log(1-a) + rnorm(n, sd = tau2)
  mat[, 2] <- qnorm(1/(1+exp(-aux2)), mean = dat$m, sd = sqrt(dat$s2))
  mat

}

