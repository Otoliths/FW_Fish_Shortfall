library(brms)

# ------------------------------------------------------------
# Stable Gompertz helpers in R
# hazard: h(t) = mu * exp(gamma * t)
# ------------------------------------------------------------

gompertz_lpdf <- function(t, mu, gamma, eps = 1e-8) {
  out <- numeric(length(t))
  near0 <- abs(gamma) < eps
  
  # exponential limit when gamma -> 0
  if (any(near0)) {
    out[near0] <- log(mu[near0]) - mu[near0] * t[near0]
  }
  
  if (any(!near0)) {
    g  <- gamma[!near0]
    m  <- mu[!near0]
    tt <- t[!near0]
    out[!near0] <- log(m) + g * tt - (m / g) * (exp(g * tt) - 1)
  }
  
  out
}

gompertz_lccdf <- function(t, mu, gamma, eps = 1e-8) {
  out <- numeric(length(t))
  near0 <- abs(gamma) < eps
  
  # exponential limit
  if (any(near0)) {
    out[near0] <- -mu[near0] * t[near0]
  }
  
  if (any(!near0)) {
    g  <- gamma[!near0]
    m  <- mu[!near0]
    tt <- t[!near0]
    out[!near0] <- -(m / g) * (exp(g * tt) - 1)
  }
  
  out
}

gompertz_surv <- function(t, mu, gamma, eps = 1e-8) {
  exp(gompertz_lccdf(t, mu, gamma, eps = eps))
}

gompertz_rng <- function(mu, gamma, eps = 1e-8) {
  n <- length(mu)
  u <- runif(n)
  out <- numeric(n)
  
  near0 <- abs(gamma) < eps
  
  # exponential limit
  if (any(near0)) {
    out[near0] <- -log(1 - u[near0]) / mu[near0]
  }
  
  if (any(!near0)) {
    g <- gamma[!near0]
    m <- mu[!near0]
    uu <- u[!near0]
    
    val <- 1 + (-g / m) * log(1 - uu)
    # if numerical issue occurs, set to Inf
    tmp <- rep(Inf, length(val))
    ok <- val > 0
    tmp[ok] <- log(val[ok]) / g[ok]
    out[!near0] <- tmp
  }
  
  out
}

# E[T] = integral_0^Inf S(t) dt
gompertz_mean <- function(mu, gamma, eps = 1e-8, upper = 500) {
  n <- length(mu)
  out <- numeric(n)
  
  near0 <- abs(gamma) < eps
  
  # exponential limit
  if (any(near0)) {
    out[near0] <- 1 / mu[near0]
  }
  
  if (any(!near0)) {
    idx <- which(!near0)
    out[idx] <- vapply(idx, function(j) {
      integrate(
        f = function(tt) gompertz_surv(tt, mu[j], gamma[j], eps = eps),
        lower = 0,
        upper = upper,
        subdivisions = 1000L,
        rel.tol = 1e-7
      )$value
    }, numeric(1))
  }
  
  out
}


# Stan-side functions
stan_funs <- "
real gompertz_lpdf(real t, real mu, real gamma) {
  if (abs(gamma) < 1e-8) {
    // exponential limit
    return log(mu) - mu * t;
  } else {
    return log(mu) + gamma * t - (mu / gamma) * (exp(gamma * t) - 1);
  }
}

real gompertz_lcdf(real t, real mu, real gamma) {
  if (abs(gamma) < 1e-8) {
    return log1m_exp(-mu * t);
  } else {
    return log1m_exp(-(mu / gamma) * (exp(gamma * t) - 1));
  }
}

real gompertz_lccdf(real t, real mu, real gamma) {
  if (abs(gamma) < 1e-8) {
    return -mu * t;
  } else {
    return -(mu / gamma) * (exp(gamma * t) - 1);
  }
}

real gompertz_rng(real mu, real gamma) {
  real u = uniform_rng(0, 1);
  real val;

  if (abs(gamma) < 1e-8) {
    return -log1m(u) / mu;
  } else {
    val = 1 - (gamma / mu) * log1m(u);
    if (val <= 0) return positive_infinity();
    return log(val) / gamma;
  }
}
"

gompertz_family <- custom_family(
  "gompertz",
  dpars = c("mu", "gamma"),
  links = c("log", "identity"),
  lb = c(0, NA),
  type = "real",
  
  log_lik = function(i, prep) {
    mu    <- brms::get_dpar(prep, "mu", i = i)
    gamma <- brms::get_dpar(prep, "gamma", i = i)
    t     <- prep$data$Y[i]
    
    # default: uncensored observation
    cens <- 0L
    
    if (!is.null(prep$data$cens) && length(prep$data$cens) >= i) {
      cens <- prep$data$cens[i]
    } else if (!is.null(prep$cens) && length(prep$cens) >= i) {
      cens <- prep$cens[i]
    }
    
    if (cens == 0) {
      gompertz_lpdf(t, mu, gamma)
    } else {
      gompertz_lccdf(t, mu, gamma)
    }
  },
  
  posterior_predict = function(i, prep, ...) {
    mu    <- brms::get_dpar(prep, "mu", i = i)
    gamma <- brms::get_dpar(prep, "gamma", i = i)
    gompertz_rng(mu, gamma)
  },
  
  posterior_epred = function(prep) {
    mu    <- brms::get_dpar(prep, "mu")
    gamma <- brms::get_dpar(prep, "gamma")
    
    S <- nrow(mu)
    N <- ncol(mu)
    out <- matrix(NA_real_, nrow = S, ncol = N)
    
    for (s in seq_len(S)) {
      out[s, ] <- gompertz_mean(mu[s, ], gamma[s, ])
    }
    out
  }
)

gompertz_family_pos <- custom_family(
  "gompertz",
  dpars = c("mu", "gamma"),
  links = c("log", "log"),
  lb = c(0, 0),
  type = "real",
  log_lik = gompertz_family$log_lik,
  posterior_predict = gompertz_family$posterior_predict,
  posterior_epred = gompertz_family$posterior_epred
)