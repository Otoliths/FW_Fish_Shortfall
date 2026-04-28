# To maintain the hierarchical dependency among the three knowledge shortfalls, 
# all posterior probabilities were computed with respect to a common knowledge 
# horizon (year = 2025), corresponding to the maximum year available across all datasets. 
# Using a unified cutoff ensures that Darwinian and Wallacean shortfalls remain properly nested within the Linnaean shortfall

compute_Linnaean_prob <- function(
    fit,
    year          = 2025,
    data          = NULL,
    origin_year   = 1758,
    probs         = c(0.025, 0.5, 0.975),
    include_data  = TRUE,
    re_formula    = NULL,
    draws         = NULL,     # optional number of posterior draws to use
    return_draws  = FALSE     # if TRUE, also return full p_draws matrix
) {
  # ------------------------------------------------------------
  # 0. Basic checks
  # ------------------------------------------------------------
  if (!inherits(fit, "brmsfit")) {
    stop("`fit` must be a brmsfit object.")
  }
  if (length(year) != 1L) {
    stop("`year` must be a single numeric value.")
  }
  if (length(probs) != 3L) {
    stop("`probs` must be a numeric vector of length 3, e.g. c(0.025, 0.5, 0.975).")
  }
  
  t <- year - origin_year
  if (t <= 0) {
    stop("`year` must be greater than `origin_year`.")
  }
  
  # ------------------------------------------------------------
  # 1. Prepare prediction data
  #    If `data` is NULL, use the model frame used for fitting.
  # ------------------------------------------------------------
  if (is.null(data)) {
    data <- model.frame(fit)
  }
  N <- nrow(data)
  if (N == 0L) {
    stop("`data` has zero rows; nothing to predict.")
  }
  
  # ------------------------------------------------------------
  # 2. Posterior draws of mu and shape via posterior_epred
  #    Weibull：
  #      y ~ Weibull(shape = alpha, scale = s)
  #      mu = E[y] = s * gamma(1 + 1 / alpha)
  #    
  #      s = mu / gamma(1 + 1 / shape)
  # ------------------------------------------------------------
  mu_full <- brms::posterior_epred(
    fit,
    newdata    = data,
    dpar       = "mu",
    re_formula = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
  # mu_full: n_draws_full × N
  
  shape_full <- brms::posterior_epred(
    fit,
    newdata    = data,
    dpar       = "shape",
    re_formula = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
  # shape_full: n_draws_full × N
  
  if (!all(dim(mu_full) == dim(shape_full))) {
    stop("Dimensions of mu_full and shape_full do not match.")
  }
  
  n_draws_full <- nrow(mu_full)
  if (n_draws_full < 1L) {
    stop("No posterior draws found in `fit`.")
  }
  
  # ------------------------------------------------------------
  # 3. If `draws` is specified, thin both mu and shape consistently
  # ------------------------------------------------------------
  if (!is.null(draws)) {
    if (!is.numeric(draws) || length(draws) != 1L || draws <= 0) {
      stop("`draws` must be a positive integer or NULL.")
    }
    draws <- as.integer(draws)
    
    if (draws > n_draws_full) {
      warning("Requested more draws than available; using all available draws instead.")
      draws_idx <- seq_len(n_draws_full)
    } else {
      draws_idx <- seq_len(draws)
    }
    
    mu_draws    <- mu_full[draws_idx, , drop = FALSE]
    shape_draws <- shape_full[draws_idx, , drop = FALSE]
  } else {
    mu_draws    <- mu_full
    shape_draws <- shape_full
  }
  
  n_draws <- nrow(mu_draws)
  if (!all(dim(mu_draws) == dim(shape_draws))) {
    stop("Dimensions of mu_draws and shape_draws do not match after thinning.")
  }
  
  # ------------------------------------------------------------
  # 4. Compute survival probabilities for Weibull
  #    brms Weibull 
  #      mu = mean
  #      shape = alpha
  #      scale = mu / gamma(1 + 1 / shape)
  #
  #    
  #      S(t) = P(T > t) = exp(-(t / scale)^shape)
  # ------------------------------------------------------------
  if (any(!is.finite(mu_draws)) || any(!is.finite(shape_draws))) {
    stop("Non-finite values detected in posterior draws of mu or shape.")
  }
  if (any(mu_draws <= 0, na.rm = TRUE)) {
    stop("All posterior draws of `mu` must be > 0 for Weibull.")
  }
  if (any(shape_draws <= 0, na.rm = TRUE)) {
    stop("All posterior draws of `shape` must be > 0 for Weibull.")
  }
  
  scale_draws <- mu_draws / gamma(1 + 1 / shape_draws)
  
  p_draws <- exp(- (t / scale_draws)^shape_draws)
  # p_draws: n_draws × N
  
  # ------------------------------------------------------------
  # 5. Summarise across draws for each observation
  # ------------------------------------------------------------
  prob_mean   <- apply(p_draws, 2, mean)
  prob_median <- apply(p_draws, 2, median)
  prob_lower  <- apply(p_draws, 2, quantile, probs = probs[1])
  prob_upper  <- apply(p_draws, 2, quantile, probs = probs[3])
  
  summary_df <- data.frame(
    prob_undesc_mean   = prob_mean,
    prob_undesc_median = prob_median,
    prob_undesc_lower  = prob_lower,
    prob_undesc_upper  = prob_upper,
    year               = year,
    time_since_origin  = t,
    stringsAsFactors   = FALSE
  )
  
  # ------------------------------------------------------------
  # 6. Bind original data if requested
  # ------------------------------------------------------------
  if (include_data) {
    out_df <- cbind(data, summary_df)
  } else {
    out_df <- summary_df
  }
  
  # ------------------------------------------------------------
  # 7. Optionally return full p_draws for Bayesian aggregation
  # ------------------------------------------------------------
  if (return_draws) {
    p_draws <- as.matrix(p_draws)
    colnames(p_draws) <- if (!is.null(rownames(data))) rownames(data) else seq_len(N)
    
    return(list(
      summary   = out_df,
      p_draws   = p_draws,
      n_draws   = n_draws,
      year      = year,
      origin    = origin_year
    ))
  } else {
    return(out_df)
  }
}

compute_Linnaean_prob_country_gompertz <- function(
    fit,
    year          = 2025,
    data          = NULL,
    origin_year   = 1758,
    probs         = c(0.025, 0.5, 0.975),
    include_data  = TRUE,
    re_formula    = NULL,
    draws         = NULL,
    return_draws  = FALSE,
    gamma_tol     = 1e-8
) {
  
  if (!inherits(fit, "brmsfit")) stop("`fit` must be a brmsfit object.")
  if (length(year) != 1L) stop("`year` must be a single numeric value.")
  if (length(probs) != 3L) stop("`probs` must be length 3.")
  
  t <- year - origin_year
  if (t <= 0) stop("`year` must be greater than `origin_year`.")
  
  if (is.null(data)) data <- model.frame(fit)
  N <- nrow(data)
  if (N == 0L) stop("`data` has zero rows.")
  
  # posterior draws for mu
  mu_full <- brms::posterior_epred(
    fit,
    newdata           = data,
    dpar              = "mu",
    re_formula        = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
  
  # posterior draws for gamma
  gamma_full <- brms::posterior_epred(
    fit,
    newdata           = data,
    dpar              = "gamma",
    re_formula        = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
  
  if (!all(dim(mu_full) == dim(gamma_full))) {
    stop("Dimensions of mu_full and gamma_full do not match.")
  }
  
  n_draws_full <- nrow(mu_full)
  if (n_draws_full < 1L) stop("No posterior draws found.")
  
  # optional thinning/subsetting
  if (!is.null(draws)) {
    if (!is.numeric(draws) || length(draws) != 1L || draws <= 0) {
      stop("`draws` must be a positive integer or NULL.")
    }
    draws <- as.integer(draws)
    
    if (draws > n_draws_full) {
      warning("Requested more draws than available; using all draws.")
      draws_idx <- seq_len(n_draws_full)
    } else {
      draws_idx <- seq_len(draws)
    }
    
    mu_draws    <- mu_full[draws_idx, , drop = FALSE]
    gamma_draws <- gamma_full[draws_idx, , drop = FALSE]
  } else {
    mu_draws    <- mu_full
    gamma_draws <- gamma_full
  }
  
  n_draws <- nrow(mu_draws)
  
  # Gompertz survival function:
  # S(t) = exp(-(mu/gamma) * (exp(gamma*t) - 1))
  # if gamma ~ 0, use exponential limit: exp(-mu * t)
  p_draws <- matrix(NA_real_, nrow = n_draws, ncol = N)
  
  for (s in seq_len(n_draws)) {
    mu_s    <- mu_draws[s, ]
    gamma_s <- gamma_draws[s, ]
    
    near_zero <- abs(gamma_s) < gamma_tol
    
    # gamma approximately zero -> exponential limit
    if (any(near_zero)) {
      p_draws[s, near_zero] <- exp(-mu_s[near_zero] * t)
    }
    
    if (any(!near_zero)) {
      g  <- gamma_s[!near_zero]
      mu <- mu_s[!near_zero]
      p_draws[s, !near_zero] <- exp(-(mu / g) * (exp(g * t) - 1))
    }
  }
  
  summary_df <- data.frame(
    prob_undesc_mean   = apply(p_draws, 2, mean),
    prob_undesc_median = apply(p_draws, 2, median),
    prob_undesc_lower  = apply(p_draws, 2, quantile, probs = probs[1]),
    prob_undesc_upper  = apply(p_draws, 2, quantile, probs = probs[3]),
    year               = year,
    time_since_origin  = t,
    stringsAsFactors   = FALSE
  )
  
  if (include_data) {
    out_df <- dplyr::bind_cols(data, summary_df)
  } else {
    out_df <- summary_df
  }
  
  if (return_draws) {
    colnames(p_draws) <- if (!is.null(rownames(data))) rownames(data) else seq_len(N)
    return(list(
      summary = out_df,
      p_draws = p_draws,
      n_draws = n_draws,
      year    = year,
      origin  = origin_year
    ))
  } else {
    return(out_df)
  }
}

compute_Wallacean_prob <- function(
    fit,
    year,                        # target calendar year (e.g. 2025)
    data,
    year_var       = "year_description",  # column storing species description year
    probs          = c(0.025, 0.5, 0.975),
    include_data   = TRUE,
    re_formula     = NA,         # default: marginalize over random effects
    draws          = NULL,       # optional: number of posterior draws to use
    return_draws   = FALSE,      # if TRUE, also return full p_draws
    clamp_probs    = TRUE,       # whether to clamp probabilities numerically
    clamp_eps      = 1e-12       # clamping boundary
) {
  # ------------------------------------------------------------------
  # Goal:
  #   For each record i (typically species × basin), compute
  #     p_i = Pr(T_i > t_i*),
  #   where T_i is the time to first geographic record.
  #
  #   Here t_i* = year - year_description_i
  #   - If t_i* > 0: use the Gompertz survival function
  #   - If t_i* <= 0: the species had not yet been described by the
  #     target year, so it is excluded from Wallacean calculations -> NA
  #
  # Model assumption:
  #   - brms Gompertz model for time-to-event
  #   - mu and gamma are extracted on the distributional parameter scale
  #     required by the Gompertz survival function
  # ------------------------------------------------------------------
  
  # --- Basic checks --------------------------------------------------
  if (!inherits(fit, "brmsfit")) {
    stop("`fit` must be a brmsfit object.")
  }
  if (missing(year) || length(year) != 1L || !is.numeric(year)) {
    stop("You must supply a single numeric `year`.")
  }
  if (missing(data)) {
    stop("You must supply `data` containing at least the predictors and `year_description`.")
  }
  if (!year_var %in% names(data)) {
    stop("`data` must contain a column named `", year_var, "`.")
  }
  if (length(probs) != 3L) {
    stop("`probs` must be a numeric vector of length 3.")
  }
  
  N <- nrow(data)
  if (N == 0L) {
    stop("`data` has zero rows; nothing to predict.")
  }
  
  # --- 1. Compute record-specific target time ------------------------
  year_desc <- data[[year_var]]
  if (any(is.na(year_desc))) {
    stop("`year_var` (", year_var, ") contains NA values; cannot define t_i*.")
  }
  
  t_star <- year - year_desc
  idx_pos <- t_star > 0
  
  # --- 2. Posterior draws of Gompertz parameters --------------------
  mu_draws_full <- brms::posterior_epred(
    fit,
    newdata    = data,
    dpar       = "mu",
    re_formula = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
  
  gamma_draws_full <- brms::posterior_epred(
    fit,
    newdata    = data,
    dpar       = "gamma",
    re_formula = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
  
  if (!all(dim(mu_draws_full) == dim(gamma_draws_full))) {
    stop("Dimensions of `mu` and `gamma` posterior draws do not match.")
  }
  
  n_draws_full <- nrow(mu_draws_full)
  if (n_draws_full < 1L) {
    stop("No posterior draws found in `fit`.")
  }
  
  # --- 3. Optional thinning/subsampling -----------------------------
  if (!is.null(draws)) {
    if (!is.numeric(draws) || length(draws) != 1L || draws <= 0) {
      stop("`draws` must be a positive integer or NULL.")
    }
    draws <- as.integer(draws)
    
    if (draws > n_draws_full) {
      warning("Requested more draws than available; using all draws instead.")
      draws_idx <- seq_len(n_draws_full)
    } else {
      draws_idx <- sort(sample(seq_len(n_draws_full), size = draws))
    }
    
    mu_draws    <- mu_draws_full[draws_idx, , drop = FALSE]
    gamma_draws <- gamma_draws_full[draws_idx, , drop = FALSE]
  } else {
    mu_draws    <- mu_draws_full
    gamma_draws <- gamma_draws_full
  }
  
  if (!all(dim(mu_draws) == dim(gamma_draws))) {
    stop("Dimensions of mu_draws and gamma_draws do not match after thinning.")
  }
  
  n_draws <- nrow(mu_draws)
  
  # --- 4. Compute posterior survival probabilities ------------------
  p_draws <- matrix(NA_real_, nrow = n_draws, ncol = N)
  
  if (any(idx_pos)) {
    t_pos     <- t_star[idx_pos]
    mu_pos    <- mu_draws[, idx_pos, drop = FALSE]
    gamma_pos <- gamma_draws[, idx_pos, drop = FALSE]
    
    # Expand t to a draw-by-observation matrix
    t_mat <- matrix(t_pos, nrow = n_draws, ncol = length(t_pos), byrow = TRUE)
    
    # Use the validated Gompertz survival helper on vectorized inputs
    p_pos <- gompertz_surv(
      t     = c(t_mat),
      mu    = c(mu_pos),
      gamma = c(gamma_pos)
    )
    dim(p_pos) <- dim(mu_pos)
    
    if (clamp_probs) {
      p_pos <- pmin(pmax(p_pos, clamp_eps), 1 - clamp_eps)
    }
    
    p_draws[, idx_pos] <- p_pos
  }
  
  # --- 5. Safe summary for each record ------------------------------
  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x)
  }
  
  safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    stats::median(x)
  }
  
  safe_quantile <- function(x, prob) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    as.numeric(stats::quantile(x, probs = prob, names = FALSE))
  }
  
  prob_nongeoloc_mean   <- apply(p_draws, 2, safe_mean)
  prob_nongeoloc_median <- apply(p_draws, 2, safe_median)
  prob_nongeoloc_lower  <- apply(p_draws, 2, safe_quantile, prob = probs[1])
  prob_nongeoloc_upper  <- apply(p_draws, 2, safe_quantile, prob = probs[3])
  
  summary_df <- data.frame(
    prob_nongeoloc_mean   = prob_nongeoloc_mean,
    prob_nongeoloc_median = prob_nongeoloc_median,
    prob_nongeoloc_lower  = prob_nongeoloc_lower,
    prob_nongeoloc_upper  = prob_nongeoloc_upper,
    year                  = year,
    stringsAsFactors      = FALSE
  )
  
  # --- 6. Bind original data if requested ---------------------------
  if (include_data) {
    out_df <- cbind(data, summary_df)
  } else {
    out_df <- summary_df
  }
  
  # --- 7. Return full posterior draws if requested ------------------
  if (return_draws) {
    colnames(p_draws) <- if (!is.null(rownames(data))) rownames(data) else seq_len(N)
    return(list(
      summary    = out_df,
      p_draws    = p_draws,
      n_draws    = n_draws,
      year       = year,
      re_formula = re_formula
    ))
  } else {
    return(out_df)
  }
}


compute_Darwinian_prob <- function(
    fit,
    year,                        # target calendar year (e.g. 2025)
    data,
    year_var       = "year_description",  # column storing species description year
    probs          = c(0.025, 0.5, 0.975),
    include_data   = TRUE,
    re_formula     = NA,         # default: marginalize over random effects (recommended)
    draws          = NULL,       # optional: number of posterior draws to use
    return_draws   = FALSE,      # if TRUE, also return full p_draws
    clamp_probs    = TRUE,
    clamp_eps      = 1e-12
) {
  # ------------------------------------------------------------------
  # Goal:
  #   Darwinian shortfall: probability that a species still has no sequence
  #   by the target year `year`.
  #
  # Model assumption:
  #   - brms Gompertz time-to-event model:
  #       time | cens(1 - event) ~ ...
  #       family = gompertz, with links mu = log and gamma = identity
  #
  #   - time:
  #       if sequenced:   time = year_sequence - year_description
  #       if unsequenced: time = cutoff_year - year_description (right-censored)
  #
  #   Let T_i be the time from description to first sequencing for species i.
  #   For a target year `year`:
  #       t_i* = year - year_description_i
  #       p_i  = Pr(T_i > t_i*)
  #
  #   Gompertz parameterization:
  #       hazard(t) = mu_i * exp(gamma_i * t)
  #       S(t) = exp(-(mu_i / gamma_i) * (exp(gamma_i * t) - 1)), gamma_i != 0
  #       S(t) = exp(-mu_i * t),                                  gamma_i -> 0
  #
  #   Note:
  #     If t_i* <= 0 (the target year is earlier than or equal to the
  #     description year), the species had not yet been described by `year`,
  #     so it should not be counted in the Darwinian shortfall -> return NA.
  # ------------------------------------------------------------------
  
  # --- Basic checks --------------------------------------------------
  if (!inherits(fit, "brmsfit")) {
    stop("`fit` must be a brmsfit object.")
  }
  if (missing(year) || length(year) != 1L || !is.numeric(year)) {
    stop("You must supply a single numeric `year`.")
  }
  if (missing(data)) {
    stop("You must supply `data` containing at least the predictors and `year_description`.")
  }
  if (!year_var %in% names(data)) {
    stop("`data` must contain a column named `", year_var, "`.")
  }
  
  N <- nrow(data)
  if (N == 0L) {
    stop("`data` has zero rows; nothing to predict.")
  }
  if (length(probs) != 3L) {
    stop("`probs` must be a numeric vector of length 3, e.g. c(0.025, 0.5, 0.975).")
  }
  
  # --- 1. Compute species-specific t_star ----------------------------
  year_desc <- data[[year_var]]
  if (any(is.na(year_desc))) {
    stop("`year_var` (", year_var, ") contains NA values; cannot define t_i*.")
  }
  
  t_star <- year - year_desc           # may be <= 0
  idx_pos <- t_star > 0                # only these species had been described by `year`
  idx_nonpos <- !idx_pos               # not yet described by `year` -> excluded from Darwinian shortfall
  
  # --- 2. Obtain posterior draws of mu -------------------------------
  # dpar = "mu" returns the Gompertz mu parameter on the response scale
  mu_draws_full <- brms::posterior_epred(
    fit,
    newdata    = data,
    dpar       = "mu",
    re_formula = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )  # dimensions: n_draws_full × N
  
  n_draws_full <- nrow(mu_draws_full)
  if (n_draws_full < 1L) {
    stop("No posterior draws found in `fit`.")
  }
  
  # --- 3. Obtain posterior draws of gamma ----------------------------
  gamma_draws_full <- brms::posterior_epred(
    fit,
    newdata    = data,
    dpar       = "gamma",
    re_formula = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )  # dimensions: n_draws_full × N
  
  if (!all(dim(mu_draws_full) == dim(gamma_draws_full))) {
    stop("Dimensions of `mu` and `gamma` posterior draws do not match.")
  }
  
  # --- 4. Optionally restrict the number of posterior draws ---------
  if (!is.null(draws)) {
    if (!is.numeric(draws) || length(draws) != 1L || draws <= 0) {
      stop("`draws` must be a positive integer or NULL.")
    }
    draws <- as.integer(draws)
    
    if (draws > n_draws_full) {
      warning("Requested more draws than available; using all draws instead.")
      draws_idx <- seq_len(n_draws_full)
    } else {
      # Use a random subset rather than the first `draws` draws
      draws_idx <- sort(sample(seq_len(n_draws_full), size = draws))
    }
    
    mu_draws    <- mu_draws_full[draws_idx, , drop = FALSE]
    gamma_draws <- gamma_draws_full[draws_idx, , drop = FALSE]
  } else {
    mu_draws    <- mu_draws_full
    gamma_draws <- gamma_draws_full
  }
  
  n_draws <- nrow(mu_draws)
  if (!all(dim(mu_draws) == dim(gamma_draws))) {
    stop("Dimensions of mu_draws and gamma_draws do not match after thinning.")
  }
  
  # --- 5. Compute p_i = Pr(T_i > t_i*) -------------------------------
  # Initialize everything as NA; records with t_star <= 0 remain NA
  p_draws <- matrix(NA_real_, nrow = n_draws, ncol = N)
  
  if (any(idx_pos)) {
    t_pos     <- t_star[idx_pos]                           # > 0
    mu_pos    <- mu_draws[, idx_pos, drop = FALSE]        # n_draws × n_pos
    gamma_pos <- gamma_draws[, idx_pos, drop = FALSE]     # n_draws × n_pos
    
    # For simplicity and to stay close to the original structure,
    # compute survival draw by draw.
    for (m in seq_len(n_draws)) {
      mu_m    <- mu_pos[m, ]
      gamma_m <- gamma_pos[m, ]
      
      p_m <- gompertz_surv(
        t     = t_pos,
        mu    = mu_m,
        gamma = gamma_m
      )
      
      if (clamp_probs) {
        p_m <- pmin(pmax(p_m, clamp_eps), 1 - clamp_eps)
      }
      
      p_draws[m, idx_pos] <- p_m
    }
  }
  
  # --- 6. Summarize posterior distributions for each record ---------
  prob_noseq_mean   <- apply(p_draws, 2, mean,   na.rm = TRUE)
  prob_noseq_median <- apply(p_draws, 2, median, na.rm = TRUE)
  prob_noseq_lower  <- apply(p_draws, 2, quantile, probs = probs[1], na.rm = TRUE)
  prob_noseq_upper  <- apply(p_draws, 2, quantile, probs = probs[3], na.rm = TRUE)
  
  summary_df <- data.frame(
    prob_noseq_mean   = prob_noseq_mean,
    prob_noseq_median = prob_noseq_median,
    prob_noseq_lower  = prob_noseq_lower,
    prob_noseq_upper  = prob_noseq_upper,
    year              = year,
    stringsAsFactors  = FALSE
  )
  
  # --- 7. Bind original data if requested ---------------------------
  if (include_data) {
    out_df <- cbind(data, summary_df)
  } else {
    out_df <- summary_df
  }
  
  # --- 8. Return full p_draws matrix if requested -------------------
  if (return_draws) {
    colnames(p_draws) <- if (!is.null(rownames(data))) rownames(data) else seq_len(N)
    return(list(
      summary    = out_df,
      p_draws    = p_draws,
      n_draws    = n_draws,
      year       = year,
      re_formula = re_formula
    ))
  } else {
    return(out_df)
  }
}


