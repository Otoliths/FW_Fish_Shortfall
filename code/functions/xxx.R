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
  
   mu_full <- brms::posterior_epred(
    fit,
    newdata    = data,
    dpar       = "mu",
    re_formula = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
  # mu_full: n_draws_full × N
  
  sigma_full <- brms::posterior_epred(
    fit,
    newdata    = data,
    dpar       = "sigma",
    re_formula = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
 
  
  if (!all(dim(mu_full) == dim(sigma_full))) {
    stop("Dimensions of mu_full and sigma_full do not match.")
  }
  
  n_draws_full <- nrow(mu_full)
  if (n_draws_full < 1L) {
    stop("No posterior draws found in `fit`.")
  }
  
  # ------------------------------------------------------------
  # 3. If `draws` is specified, thin both mu and sigma consistently
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
    sigma_draws <- sigma_full[draws_idx, , drop = FALSE]
  } else {
    mu_draws    <- mu_full
    sigma_draws <- sigma_full
  }
  
  n_draws <- nrow(mu_draws)
  if (!all(dim(mu_draws) == dim(sigma_draws))) {
    stop("Dimensions of mu_draws and sigma_draws do not match after thinning.")
  }
  
  
  p_draws <- plnorm(
    q          = t,
    meanlog    = mu_draws,
    sdlog      = sigma_draws,
    lower.tail = FALSE
  )
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

compute_Linnaean_prob_country <- function(
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
  if (!inherits(fit, "brmsfit")) stop("`fit` must be a brmsfit object.")
  if (length(year) != 1L) stop("`year` must be a single numeric value.")
  if (length(probs) != 3L) stop("`probs` must be a numeric vector of length 3, e.g. c(0.025, 0.5, 0.975).")
  
  t <- year - origin_year
  if (t <= 0) stop("`year` must be greater than `origin_year`.")
  
  # ------------------------------------------------------------
  # 1. Prepare prediction data
  # ------------------------------------------------------------
  if (is.null(data)) data <- model.frame(fit)
  N <- nrow(data)
  if (N == 0L) stop("`data` has zero rows; nothing to predict.")
  
  # ------------------------------------------------------------
  # 2. Posterior draws of mu and shape via posterior_epred (Weibull)
  #    T ~ Weibull(shape, scale = mu)
  # ------------------------------------------------------------
  mu_full <- brms::posterior_epred(
    fit,
    newdata           = data,
    dpar              = "mu",
    re_formula        = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
  
  shape_full <- brms::posterior_epred(
    fit,
    newdata           = data,
    dpar              = "shape",
    re_formula        = re_formula,
    allow_new_levels  = TRUE,
    sample_new_levels = "gaussian"
  )
  
  if (!all(dim(mu_full) == dim(shape_full))) stop("Dimensions of mu_full and shape_full do not match.")
  
  n_draws_full <- nrow(mu_full)
  if (n_draws_full < 1L) stop("No posterior draws found in `fit`.")
  
  # ------------------------------------------------------------
  # 3. Optional thinning of posterior draws
  # ------------------------------------------------------------
  if (!is.null(draws)) {
    if (!is.numeric(draws) || length(draws) != 1L || draws <= 0) stop("`draws` must be a positive integer or NULL.")
    draws <- as.integer(draws)
    
    if (draws > n_draws_full) {
      warning("Requested more draws than available; using all available draws instead.")
      draws_idx <- seq_len(n_draws_full)
    } else {
      draws_idx <- seq_len(draws)  # deterministic; replace with sample() if desired
    }
    
    mu_draws    <- mu_full[draws_idx, , drop = FALSE]
    shape_draws <- shape_full[draws_idx, , drop = FALSE]
  } else {
    mu_draws    <- mu_full
    shape_draws <- shape_full
  }
  
  n_draws <- nrow(mu_draws)
  if (!all(dim(mu_draws) == dim(shape_draws))) stop("Dimensions of mu_draws and shape_draws do not match after thinning.")
  
  # ------------------------------------------------------------
  # 4. Survival probabilities: S(t)=P(T>t)
  # ------------------------------------------------------------
  p_draws <- pweibull(
    q          = t,
    shape      = shape_draws,
    scale      = mu_draws,
    lower.tail = FALSE
  )
  
  # ------------------------------------------------------------
  # 5. Summarise across draws for each observation
  # ------------------------------------------------------------
  summary_df <- data.frame(
    prob_undesc_mean   = apply(p_draws, 2, mean),
    prob_undesc_median = apply(p_draws, 2, median),
    prob_undesc_lower  = apply(p_draws, 2, quantile, probs = probs[1]),
    prob_undesc_upper  = apply(p_draws, 2, quantile, probs = probs[3]),
    year               = year,
    time_since_origin  = t,
    stringsAsFactors   = FALSE
  )
  
  # ------------------------------------------------------------
  # 6. Bind original data (avoid duplicate-name issues from cbind)
  # ------------------------------------------------------------
  if (include_data) {
    out_df <- dplyr::bind_cols(data, summary_df)
  } else {
    out_df <- summary_df
  }
  
  # ------------------------------------------------------------
  # 7. Optionally return full p_draws matrix
  # ------------------------------------------------------------
  if (return_draws) {
    p_draws <- as.matrix(p_draws)
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
    year_var       = "year_description",  
    probs          = c(0.025, 0.5, 0.975),
    include_data   = TRUE,
    re_formula     = NA,         
    draws          = NULL,       
    return_draws   = FALSE,      
    clamp_probs    = TRUE,       
    clamp_eps      = 1e-12       
) {
 
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
  
  # --- 1. 计算物种/记录的目标时间 t* --------------------------------
  year_desc <- data[[year_var]]
  if (any(is.na(year_desc))) {
    stop("`year_var` (", year_var, ") contains NA values; cannot define t_i*.")
  }
  
  t_star <- year - year_desc      
  idx_pos <- t_star > 0          
  idx_nonpos <- !idx_pos          
  
  # --- 2. posterior_linpred  ----------------------
 
  mu_draws_full <- brms::posterior_linpred(
    fit,
    newdata    = data,
    transform  = FALSE,    
    re_formula = re_formula
  )  
  
  n_draws_full <- nrow(mu_draws_full)
  if (n_draws_full < 1L) {
    stop("No posterior draws found in `fit`.")
  }
  
  
  draws_df <- brms::as_draws_df(fit)
  sigma_candidates <- grep("^sigma($|_)", names(draws_df), value = TRUE)
  if (length(sigma_candidates) != 1L) {
    stop(
      "Could not uniquely identify `sigma` parameter in the brms fit. ",
      "Found: ", paste(sigma_candidates, collapse = ", ")
    )
  }
  sigma_full <- draws_df[[sigma_candidates]]
  
  if (length(sigma_full) != n_draws_full) {
    stop("Length of `sigma` draws does not match the number of mu draws.")
  }
  
 
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
    
    mu_draws  <- mu_draws_full[draws_idx, , drop = FALSE]
    sigma_vec <- sigma_full[draws_idx]
  } else {
    mu_draws  <- mu_draws_full
    sigma_vec <- sigma_full
  }
  
  n_draws <- nrow(mu_draws)
  if (n_draws != length(sigma_vec)) {
    stop("Number of posterior draws for mu and sigma do not match after thinning.")
  }
  
  
  p_draws <- matrix(NA_real_, nrow = n_draws, ncol = N)
  
  if (any(idx_pos)) {
    t_pos <- t_star[idx_pos]           
    log_t <- log(t_pos)                 
    mu_pos <- mu_draws[, idx_pos, drop = FALSE]  # n_draws × n_pos
    
    p_pos <- 1 - pnorm(z_mat)
    
    if (clamp_probs) {
      p_pos <- pmin(pmax(p_pos, clamp_eps), 1 - clamp_eps)
    }
    
    p_draws[, idx_pos] <- p_pos
  }
  
  
  prob_nongeoloc_mean   <- apply(p_draws, 2, mean,   na.rm = TRUE)
  prob_nongeoloc_median <- apply(p_draws, 2, median, na.rm = TRUE)
  prob_nongeoloc_lower  <- apply(p_draws, 2, quantile, probs = probs[1], na.rm = TRUE)
  prob_nongeoloc_upper  <- apply(p_draws, 2, quantile, probs = probs[3], na.rm = TRUE)
  
  
  
  summary_df <- data.frame(
    prob_nongeoloc_mean   = prob_nongeoloc_mean,
    prob_nongeoloc_median = prob_nongeoloc_median,
    prob_nongeoloc_lower  = prob_nongeoloc_lower,
    prob_nongeoloc_upper  = prob_nongeoloc_upper,
    year                  = year,
    stringsAsFactors      = FALSE
  )
  
  
  if (include_data) {
    out_df <- cbind(data, summary_df)
  } else {
    out_df <- summary_df
  }
  
  
  if (return_draws) {
    colnames(p_draws) <- if (!is.null(rownames(data))) rownames(data) else seq_len(N)
    return(list(
      summary   = out_df,
      p_draws   = p_draws,
      n_draws   = n_draws,
      year      = year,
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
    year_var       = "year_description",  
    probs          = c(0.025, 0.5, 0.975),
    include_data   = TRUE,
    re_formula     = NA,         
    draws          = NULL,       
    return_draws   = FALSE,      
    clamp_probs    = TRUE,
    clamp_eps      = 1e-12
) {
  
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
  
  
  year_desc <- data[[year_var]]
  if (any(is.na(year_desc))) {
    stop("`year_var` (", year_var, ") contains NA values; cannot define t_i*.")
  }
  
  t_star <- year - year_desc           
  idx_pos <- t_star > 0               
  idx_nonpos <- !idx_pos             
  
  
  mu_draws_full <- brms::posterior_linpred(
    fit,
    newdata    = data,
    transform  = TRUE,        
    re_formula = re_formula
  )  
  
  n_draws_full <- nrow(mu_draws_full)
  if (n_draws_full < 1L) {
    stop("No posterior draws found in `fit`.")
  }
  
  
  draws_df <- brms::as_draws_df(fit)
  shape_candidates <- grep("^shape($|_)", names(draws_df), value = TRUE)
  if (length(shape_candidates) != 1L) {
    stop(
      "Could not uniquely identify `shape` parameter in the brms fit. ",
      "Found: ", paste(shape_candidates, collapse = ", ")
    )
  }
  shape_full <- draws_df[[shape_candidates]]
  
  if (length(shape_full) != n_draws_full) {
    stop("Length of `shape` draws does not match the number of mu draws.")
  }
  
  
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
    
    mu_draws   <- mu_draws_full[draws_idx, , drop = FALSE]
    shape_vec  <- shape_full[draws_idx]
  } else {
    mu_draws   <- mu_draws_full
    shape_vec  <- shape_full
  }
  
  n_draws <- nrow(mu_draws)
  if (n_draws != length(shape_vec)) {
    stop("Number of posterior draws for mu and shape do not match after thinning.")
  }
  
  
  p_draws <- matrix(NA_real_, nrow = n_draws, ncol = N)
  
  if (any(idx_pos)) {
    t_pos   <- t_star[idx_pos]          # > 0
    mu_pos  <- mu_draws[, idx_pos, drop = FALSE]  # n_draws × n_pos
    
   
    for (m in seq_len(n_draws)) {
      shape_m  <- shape_vec[m]
      mu_m     <- mu_pos[m, ]
      scale_m  <- mu_m / shape_m
      
      p_m      <- 1 - pgamma(t_pos, shape = shape_m, scale = scale_m)
      
      if (clamp_probs) {
        p_m <- pmin(pmax(p_m, clamp_eps), 1 - clamp_eps)
      }
      
      p_draws[m, idx_pos] <- p_m
    }
  }
  
  
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
  
 
  if (include_data) {
    out_df <- cbind(data, summary_df)
  } else {
    out_df <- summary_df
  }
  
  
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

compute_Darwinian_prob_country <- function(
    fit,
    year,                        # target calendar year (e.g. 2025)
    data,
    year_var       = "year_description",  # column: species' description year
    probs          = c(0.025, 0.5, 0.975),
    include_data   = TRUE,
    re_formula     = NA,         # default: marginalize over random effects
    draws          = NULL,       # optional: number of posterior draws to use
    return_draws   = FALSE,      # TRUE = also return full p_draws matrix
    clamp_probs    = TRUE,
    clamp_eps      = 1e-12
) {
  # ------------------------------------------------------------------
  # Goal:
  #   Darwinian shortfall at species level:
  #   Probability that a species still lacks any molecular sequence
  #   by target year `year`.
  #
  # Model assumption:
  #   - brms Weibull AFT model:
  #       time | cens(1 - event) ~ ...
  #       family = weibull()
  #
  #   - time:
  #       if sequenced:   time = year_sequence - year_description
  #       if unsequenced: time = cutoff_year   - year_description  (right-censored)
  #
  #   Let T_i be "time from description to first sequence".
  #   For target calendar year `year`:
  #       t_i* = year - year_description_i
  #       p_i  = Pr(T_i > t_i*)
  #
  #   Weibull parameterisation in brms:
  #       T_i ~ Weibull(shape = alpha, scale = s_i)
  #       with mean  mu_i  and shape alpha:
  #           s_i = mu_i / Gamma(1 + 1 / alpha)
  #
  #   Note:
  #     If t_i* <= 0 (target year <= description year), the species is not
  #     yet formally described at `year` → it does not contribute to the
  #     Darwinian shortfall by definition; we return NA for these records.
  # ------------------------------------------------------------------
  
  # --- Basic checks --------------------------------------------------
  if (!inherits(fit, "brmsfit")) {
    stop("`fit` must be a brmsfit object.")
  }
  if (missing(year) || length(year) != 1L || !is.numeric(year)) {
    stop("You must supply a single numeric `year`.")
  }
  if (missing(data)) {
    stop("You must supply `data` containing predictors and `year_description`.")
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
  
  # --- 1. Species-specific time horizon t_star -----------------------
  year_desc <- data[[year_var]]
  if (any(is.na(year_desc))) {
    stop("`year_var` (", year_var, ") contains NA values; cannot define t_i*.")
  }
  
  t_star    <- year - year_desc     # may be <= 0
  idx_pos   <- t_star > 0           # described before target year
  idx_nonpos <- !idx_pos            # not yet described at `year` → NA
  
  # --- 2. Posterior draws of mu (Weibull mean) -----------------------
  # posterior_linpred(..., transform = TRUE, dpar = "mu")
  # returns draws of mu on the response scale.
  mu_draws_full <- brms::posterior_linpred(
    fit,
    newdata    = data,
    dpar       = "mu",
    transform  = TRUE,       # inverse link, gives mu on time scale
    re_formula = re_formula
  )  # n_draws_full × N
  
  n_draws_full <- nrow(mu_draws_full)
  if (n_draws_full < 1L) {
    stop("No posterior draws found in `fit`.")
  }
  
  # --- 3. Posterior draws of shape parameter -------------------------
  # For Weibull in brms, there is a global "shape" parameter with its own prior.
  # We extract it from the draws data frame.
  draws_df <- brms::as_draws_df(fit)
  shape_candidates <- grep("^shape($|_)", names(draws_df), value = TRUE)
  if (length(shape_candidates) != 1L) {
    stop(
      "Could not uniquely identify `shape` parameter in the brms fit. ",
      "Found: ", paste(shape_candidates, collapse = ", ")
    )
  }
  shape_full <- draws_df[[shape_candidates]]
  
  if (length(shape_full) != n_draws_full) {
    stop("Length of `shape` draws does not match the number of mu draws.")
  }
  
  # --- 4. Optional thinning of posterior draws -----------------------
  if (!is.null(draws)) {
    if (!is.numeric(draws) || length(draws) != 1L || draws <= 0) {
      stop("`draws` must be a positive integer or NULL.")
    }
    draws <- as.integer(draws)
    
    if (draws > n_draws_full) {
      warning("Requested more draws than available; using all draws instead.")
      draws_idx <- seq_len(n_draws_full)
    } else {
      # random subset of draws for better Monte Carlo representativeness
      draws_idx <- sort(sample(seq_len(n_draws_full), size = draws))
    }
    
    mu_draws  <- mu_draws_full[draws_idx, , drop = FALSE]
    shape_vec <- shape_full[draws_idx]
  } else {
    mu_draws  <- mu_draws_full
    shape_vec <- shape_full
  }
  
  n_draws <- nrow(mu_draws)
  if (n_draws != length(shape_vec)) {
    stop("Number of posterior draws for mu and shape do not match after thinning.")
  }
  
  # --- 5. Compute p_i = Pr(T_i > t_i*) under Weibull -----------------
  # Initialise as NA; entries with t_star <= 0 remain NA.
  p_draws <- matrix(NA_real_, nrow = n_draws, ncol = N)
  
  if (any(idx_pos)) {
    t_pos  <- t_star[idx_pos]                    # > 0
    mu_pos <- mu_draws[, idx_pos, drop = FALSE]  # n_draws × n_pos
    
    # For each draw m:
    #   alpha_m  = shape_vec[m]                    (shape)
    #   mu_mi    = mu_pos[m, ]
    #   scale_mi = mu_mi / Gamma(1 + 1 / alpha_m)  (Weibull scale)
    #   p_mi     = P(T > t_pos) = pweibull(t_pos, shape = alpha_m,
    #                                      scale = scale_mi,
    #                                      lower.tail = FALSE)
    for (m in seq_len(n_draws)) {
      alpha_m <- shape_vec[m]
      mu_m    <- mu_pos[m, ]
      
      # scale s = mu / Gamma(1 + 1/alpha)
      scale_m <- mu_m / base::gamma(1 + 1 / alpha_m)
      
      p_m <- stats::pweibull(
        q          = t_pos,
        shape      = alpha_m,
        scale      = scale_m,
        lower.tail = FALSE
      )
      
      if (clamp_probs) {
        p_m <- pmin(pmax(p_m, clamp_eps), 1 - clamp_eps)
      }
      
      p_draws[m, idx_pos] <- p_m
    }
  }
  
  # --- 6. Summarise posterior probabilities per species --------------
  prob_noseq_mean   <- apply(p_draws, 2, mean,   na.rm = TRUE)
  prob_noseq_median <- apply(p_draws, 2, median, na.rm = TRUE)
  prob_noseq_lower  <- apply(p_draws, 2, quantile,
                             probs = probs[1], na.rm = TRUE)
  prob_noseq_upper  <- apply(p_draws, 2, quantile,
                             probs = probs[3], na.rm = TRUE)
  
  summary_df <- data.frame(
    prob_noseq_mean   = prob_noseq_mean,
    prob_noseq_median = prob_noseq_median,
    prob_noseq_lower  = prob_noseq_lower,
    prob_noseq_upper  = prob_noseq_upper,
    year              = year,
    stringsAsFactors  = FALSE
  )
  
  # --- 7. Optionally bind original data ------------------------------
  if (include_data) {
    out_df <- cbind(data, summary_df)
  } else {
    out_df <- summary_df
  }
  
  # --- 8. Optionally return full p_draws matrix ----------------------
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



energy_test_best_worst <- function(data,
                                   vars_test,
                                   strategy_col = "strategy",
                                   best_label  = "best",
                                   worst_label = "worst",
                                   scale_vars  = TRUE,
                                   R = 999,
                                   seed = 1,
                                   min_n = 2) {
  
  
  dat_bw <- data %>%
    dplyr::filter(.data[[strategy_col]] %in% c(best_label, worst_label)) %>%
    dplyr::select(dplyr::all_of(vars_test), dplyr::all_of(strategy_col)) %>%
    tidyr::drop_na()
  
  n_best  <- sum(dat_bw[[strategy_col]] == best_label)
  n_worst <- sum(dat_bw[[strategy_col]] == worst_label)
  
  if (n_best < min_n || n_worst < min_n) {
    return(tibble::tibble(
      n_best    = n_best,
      n_worst   = n_worst,
      statistic = NA_real_,
      p_value   = NA_real_,
      note      = "Insufficient sample size"
    ))
  }
  
 
  X <- dat_bw %>%
    dplyr::select(dplyr::all_of(vars_test)) %>%
    {
      if (scale_vars) {
        dplyr::mutate(., dplyr::across(everything(), ~ as.numeric(scale(.x))))
      } else .
    } %>%
    as.matrix()
  
  # ---- 3. Energy test ----
  set.seed(seed)
  res <- energy::eqdist.etest(
    X,
    sizes = c(n_best, n_worst),
    R = R
  )
  
  
  tibble::tibble(
    n_best    = n_best,
    n_worst   = n_worst,
    statistic = unname(res$statistic),
    p_value   = unname(res$p.value)
  )
}



