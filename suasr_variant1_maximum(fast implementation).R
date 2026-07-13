# SuASR Variant I: exact maximum infectious-burden distribution
#
# This single-file script implements the recurrent-infection formulation of
# the Markovian SuASR model with stochastic external factor states. It computes
# the exact distribution of
#
#   M = max_t { I(t) + A(t) }
#
# until epidemic extinction and validates the result with direct Gillespie
# simulation. The numerical example is Manuscript validation case 1.
#
# Required package: Matrix

.require_matrix <- function() {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package 'Matrix' is required. Install it with install.packages('Matrix').",
      call. = FALSE
    )
  }
}

.assert_integerish <- function(x, name, lower = 0L) {
  if (length(x) != 1L || !is.finite(x) || abs(x - round(x)) > 1e-10) {
    stop(sprintf("'%s' must be a finite integer.", name), call. = FALSE)
  }

  x <- as.integer(round(x))
  if (x < lower) {
    stop(sprintf("'%s' must be at least %d.", name, lower), call. = FALSE)
  }

  x
}

.expand_parameter <- function(x, n_efs, name) {
  if (length(x) == 1L) {
    x <- rep(as.numeric(x), n_efs)
  } else if (length(x) == n_efs) {
    x <- as.numeric(x)
  } else {
    stop(
      sprintf("'%s' must have length 1 or length %d.", name, n_efs),
      call. = FALSE
    )
  }

  if (any(!is.finite(x))) {
    stop(sprintf("'%s' must contain only finite values.", name), call. = FALSE)
  }

  x
}

# Standardize the generator of the external factor process. Only its
# off-diagonal entries are used; each diagonal entry is reconstructed as the
# negative row sum of the corresponding off-diagonal rates.
standardize_efs_generator <- function(Q) {
  Q <- as.matrix(Q)

  if (length(dim(Q)) != 2L || nrow(Q) != ncol(Q) || nrow(Q) < 1L) {
    stop("'Q' must be a nonempty square matrix.", call. = FALSE)
  }
  if (any(!is.finite(Q))) {
    stop("'Q' must contain only finite values.", call. = FALSE)
  }

  off_diagonal <- Q
  diag(off_diagonal) <- 0

  tolerance <- 100 * .Machine$double.eps
  if (any(off_diagonal < -tolerance)) {
    stop("The off-diagonal entries of 'Q' must be nonnegative.", call. = FALSE)
  }

  off_diagonal[off_diagonal < 0] <- 0
  Q_standard <- off_diagonal
  diag(Q_standard) <- -rowSums(off_diagonal)
  Q_standard
}

new_variant1_parameters <- function(
    N,
    p,
    beta,
    gamma,
    kappa,
    mu,
    Q
) {
  N <- .assert_integerish(N, "N", lower = 1L)
  Q <- standardize_efs_generator(Q)
  n_efs <- nrow(Q)

  p <- .expand_parameter(p, n_efs, "p")
  beta <- .expand_parameter(beta, n_efs, "beta")
  gamma <- .expand_parameter(gamma, n_efs, "gamma")
  kappa <- .expand_parameter(kappa, n_efs, "kappa")
  mu <- .expand_parameter(mu, n_efs, "mu")

  if (any(p < 0 | p > 1)) {
    stop("All entries of 'p' must lie in [0, 1].", call. = FALSE)
  }

  rate_vectors <- list(beta = beta, gamma = gamma, kappa = kappa, mu = mu)
  invalid <- vapply(rate_vectors, function(x) any(x < 0), logical(1))
  if (any(invalid)) {
    stop(
      sprintf(
        "The following rate vectors contain negative values: %s.",
        paste(names(rate_vectors)[invalid], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  structure(
    list(
      N = N,
      n_efs = n_efs,
      p = p,
      beta = beta,
      gamma = gamma,
      kappa = kappa,
      mu = mu,
      Q = Q
    ),
    class = "suasr_variant1_parameters"
  )
}

.validate_initial_state <- function(init, parameters) {
  if (is.null(names(init))) {
    stop(
      "'init' must be a named vector containing I, A, R, and F; S is optional.",
      call. = FALSE
    )
  }

  required <- c("I", "A", "R", "F")
  if (!all(required %in% names(init))) {
    stop(
      "'init' must contain named entries I, A, R, and F; S is optional.",
      call. = FALSE
    )
  }

  values <- as.numeric(init[required])
  if (any(!is.finite(values)) || any(abs(values - round(values)) > 1e-10)) {
    stop("I, A, R, and F in 'init' must be finite integers.", call. = FALSE)
  }

  values <- as.integer(round(values))
  names(values) <- required

  I <- unname(values["I"])
  A <- unname(values["A"])
  R <- unname(values["R"])
  F <- unname(values["F"])
  S <- parameters$N - I - A - R

  if (any(c(S, I, A, R) < 0L)) {
    stop("The initial compartment counts must be nonnegative and sum to N.", call. = FALSE)
  }
  if (I + A == 0L) {
    stop("The initial state must contain at least one infectious individual.", call. = FALSE)
  }
  if (F < 1L || F > parameters$n_efs) {
    stop("The initial external factor state F is outside the valid range.", call. = FALSE)
  }

  if ("S" %in% names(init)) {
    supplied_S <- as.numeric(init["S"])
    if (!is.finite(supplied_S) || abs(supplied_S - round(supplied_S)) > 1e-10) {
      stop("S in 'init' must be a finite integer.", call. = FALSE)
    }
    if (as.integer(round(supplied_S)) != S) {
      stop("The supplied initial state does not satisfy S + I + A + R = N.", call. = FALSE)
    }
  }

  c(S = S, I = I, A = A, R = R, F = F)
}

.state_key <- function(I, A, R, F) {
  paste(I, A, R, F, sep = ":")
}

# Enumerate all transient epidemic states. A state is transient when I + A > 0.
.enumerate_transient_states <- function(parameters) {
  N <- parameters$N
  rows <- vector("list", parameters$n_efs)

  for (F in seq_len(parameters$n_efs)) {
    current <- vector("list", as.integer(round(choose(N + 3L, 3L))))
    position <- 0L

    for (I in 0:N) {
      for (A in 0:(N - I)) {
        for (R in 0:(N - I - A)) {
          if (I + A == 0L) {
            next
          }

          position <- position + 1L
          S <- N - I - A - R
          current[[position]] <- c(S = S, I = I, A = A, R = R, F = F)
        }
      }
    }

    current <- current[seq_len(position)]
    rows[[F]] <- do.call(rbind, current)
  }

  states <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  rownames(states) <- NULL
  states[] <- lapply(states, as.integer)
  states$active <- states$I + states$A
  states
}

# Construct the full transient generator T and the vector q_ext of rates from
# transient states to the disease-free set. The diagonal of T contains the
# negative total rate of every possible event, including exits from any
# threshold-restricted state space used later.
build_variant1_transient_generator <- function(parameters, check = TRUE) {
  .require_matrix()

  if (!inherits(parameters, "suasr_variant1_parameters")) {
    stop(
      "'parameters' must be created by new_variant1_parameters().",
      call. = FALSE
    )
  }

  states <- .enumerate_transient_states(parameters)
  n_states <- nrow(states)

  keys <- .state_key(states$I, states$A, states$R, states$F)
  state_lookup <- setNames(seq_len(n_states), keys)

  lookup_state <- function(I, A, R, F) {
    id <- unname(state_lookup[.state_key(I, A, R, F)])
    if (length(id) != 1L || is.na(id)) {
      stop("Internal error: a destination state was not found.", call. = FALSE)
    }
    as.integer(id)
  }

  row_list <- vector("list", n_states)
  col_list <- vector("list", n_states)
  value_list <- vector("list", n_states)
  q_ext <- numeric(n_states)

  add_transition <- function(row, destination, rate, local_rows, local_cols, local_values) {
    if (rate <= 0) {
      return(list(rows = local_rows, cols = local_cols, values = local_values))
    }

    destination_I <- destination[1L]
    destination_A <- destination[2L]
    destination_R <- destination[3L]
    destination_F <- destination[4L]

    if (destination_I + destination_A == 0L) {
      q_ext[row] <<- q_ext[row] + rate
    } else {
      destination_id <- lookup_state(
        destination_I,
        destination_A,
        destination_R,
        destination_F
      )
      local_rows <- c(local_rows, row)
      local_cols <- c(local_cols, destination_id)
      local_values <- c(local_values, rate)
    }

    list(rows = local_rows, cols = local_cols, values = local_values)
  }

  for (row in seq_len(n_states)) {
    S <- states$S[row]
    I <- states$I[row]
    A <- states$A[row]
    R <- states$R[row]
    F <- states$F[row]
    infectious <- I + A

    local_rows <- integer(0)
    local_cols <- integer(0)
    local_values <- numeric(0)
    total_rate <- 0

    rate_sym_infection <- parameters$p[F] * parameters$beta[F] * S * infectious / parameters$N
    rate_asym_infection <- (1 - parameters$p[F]) * parameters$beta[F] * S * infectious / parameters$N
    rate_sym_recovery <- parameters$gamma[F] * I
    rate_asym_recovery <- parameters$kappa[F] * A
    rate_sym_removal <- parameters$mu[F] * I

    if (S > 0L && rate_sym_infection > 0) {
      update <- add_transition(
        row,
        c(I + 1L, A, R, F),
        rate_sym_infection,
        local_rows,
        local_cols,
        local_values
      )
      local_rows <- update$rows
      local_cols <- update$cols
      local_values <- update$values
      total_rate <- total_rate + rate_sym_infection
    }

    if (S > 0L && rate_asym_infection > 0) {
      update <- add_transition(
        row,
        c(I, A + 1L, R, F),
        rate_asym_infection,
        local_rows,
        local_cols,
        local_values
      )
      local_rows <- update$rows
      local_cols <- update$cols
      local_values <- update$values
      total_rate <- total_rate + rate_asym_infection
    }

    if (I > 0L && rate_sym_recovery > 0) {
      update <- add_transition(
        row,
        c(I - 1L, A, R, F),
        rate_sym_recovery,
        local_rows,
        local_cols,
        local_values
      )
      local_rows <- update$rows
      local_cols <- update$cols
      local_values <- update$values
      total_rate <- total_rate + rate_sym_recovery
    }

    if (A > 0L && rate_asym_recovery > 0) {
      update <- add_transition(
        row,
        c(I, A - 1L, R, F),
        rate_asym_recovery,
        local_rows,
        local_cols,
        local_values
      )
      local_rows <- update$rows
      local_cols <- update$cols
      local_values <- update$values
      total_rate <- total_rate + rate_asym_recovery
    }

    if (I > 0L && rate_sym_removal > 0) {
      update <- add_transition(
        row,
        c(I - 1L, A, R + 1L, F),
        rate_sym_removal,
        local_rows,
        local_cols,
        local_values
      )
      local_rows <- update$rows
      local_cols <- update$cols
      local_values <- update$values
      total_rate <- total_rate + rate_sym_removal
    }

    efs_targets <- setdiff(seq_len(parameters$n_efs), F)
    for (target_F in efs_targets) {
      rate_efs <- parameters$Q[F, target_F]
      if (rate_efs <= 0) {
        next
      }

      update <- add_transition(
        row,
        c(I, A, R, target_F),
        rate_efs,
        local_rows,
        local_cols,
        local_values
      )
      local_rows <- update$rows
      local_cols <- update$cols
      local_values <- update$values
      total_rate <- total_rate + rate_efs
    }

    local_rows <- c(local_rows, row)
    local_cols <- c(local_cols, row)
    local_values <- c(local_values, -total_rate)

    row_list[[row]] <- local_rows
    col_list[[row]] <- local_cols
    value_list[[row]] <- local_values
  }

  T <- Matrix::sparseMatrix(
    i = unlist(row_list, use.names = FALSE),
    j = unlist(col_list, use.names = FALSE),
    x = unlist(value_list, use.names = FALSE),
    dims = c(n_states, n_states),
    giveCsparse = TRUE
  )

  if (check) {
    row_balance <- as.numeric(Matrix::rowSums(T)) + q_ext
    tolerance <- 1e-10 * max(1, max(abs(Matrix::diag(T))))

    if (max(abs(row_balance)) > tolerance) {
      stop("Generator check failed: rows do not balance to zero.", call. = FALSE)
    }
    if (any(Matrix::diag(T) >= 0)) {
      stop("Generator check failed: transient diagonal entries must be negative.", call. = FALSE)
    }
  }

  list(
    T = T,
    q_ext = q_ext,
    states = states,
    state_lookup = state_lookup
  )
}

# Compute P(M <= m) for m = 0, ..., N. For a fixed threshold m, states with
# I + A > m are excluded. Because the original diagonal rates are retained,
# transitions above m act as exits to a failure state. Therefore,
#
#   -T_m^{-1} q_ext
#
# gives the probability of epidemic extinction before the threshold is
# exceeded.
exact_maximum_distribution <- function(parameters, init, check = TRUE) {
  init <- .validate_initial_state(init, parameters)
  generator <- build_variant1_transient_generator(parameters, check = check)

  initial_key <- .state_key(init["I"], init["A"], init["R"], init["F"])
  initial_id <- unname(generator$state_lookup[initial_key])

  if (length(initial_id) != 1L || is.na(initial_id)) {
    stop("The initial transient state was not found.", call. = FALSE)
  }

  initial_active <- unname(init["I"] + init["A"])
  support <- 0:parameters$N
  cdf <- numeric(length(support))

  for (threshold in initial_active:parameters$N) {
    retained <- generator$states$active <= threshold
    retained_ids <- which(retained)
    initial_local_id <- match(initial_id, retained_ids)

    T_threshold <- generator$T[retained, retained, drop = FALSE]
    q_threshold <- generator$q_ext[retained]

    hitting_probability <- -as.numeric(Matrix::solve(T_threshold, q_threshold))
    cdf[threshold + 1L] <- hitting_probability[initial_local_id]
  }

  numerical_tolerance <- 1e-9
  cdf[abs(cdf) < numerical_tolerance] <- 0
  cdf[abs(cdf - 1) < numerical_tolerance] <- 1
  cdf <- pmin(1, pmax(0, cummax(cdf)))

  if (abs(cdf[parameters$N + 1L] - 1) > 1e-7) {
    warning(
      sprintf(
        "The final CDF value is %.10f rather than 1. Check the model parameters.",
        cdf[parameters$N + 1L]
      ),
      call. = FALSE
    )
  }

  pmf <- c(cdf[1L], diff(cdf))
  pmf[abs(pmf) < numerical_tolerance] <- 0

  if (any(pmf < -1e-7)) {
    warning("Numerically negative probabilities were detected.", call. = FALSE)
  }

  pmf <- pmax(pmf, 0)
  probability_mass <- sum(pmf)
  if (abs(probability_mass - 1) < 1e-7) {
    pmf <- pmf / probability_mass
  }

  mean_maximum <- sum(support * pmf)
  variance_maximum <- sum((support - mean_maximum)^2 * pmf)

  list(
    distribution = data.frame(
      maximum = support,
      probability = pmf,
      cumulative_probability = cdf
    ),
    summary = data.frame(
      mean = mean_maximum,
      variance = variance_maximum,
      sd = sqrt(variance_maximum),
      probability_mass = sum(pmf)
    ),
    initial_state = init,
    parameters = parameters
  )
}

simulate_variant1_gillespie <- function(
    parameters,
    init,
    max_events = 1000000L
) {
  init <- .validate_initial_state(init, parameters)
  max_events <- .assert_integerish(max_events, "max_events", lower = 1L)

  S <- unname(init["S"])
  I <- unname(init["I"])
  A <- unname(init["A"])
  R <- unname(init["R"])
  F <- unname(init["F"])

  time <- 0
  maximum_infectious <- I + A
  event_count <- 0L

  while (I + A > 0L && event_count < max_events) {
    event_count <- event_count + 1L
    infectious <- I + A

    rate_asym_infection <-
      (1 - parameters$p[F]) * parameters$beta[F] * S * infectious / parameters$N
    rate_sym_infection <-
      parameters$p[F] * parameters$beta[F] * S * infectious / parameters$N
    rate_sym_recovery <- parameters$gamma[F] * I
    rate_asym_recovery <- parameters$kappa[F] * A
    rate_sym_removal <- parameters$mu[F] * I

    efs_targets <- setdiff(seq_len(parameters$n_efs), F)
    efs_rates <- parameters$Q[F, efs_targets]

    rates <- c(
      asym_infection = rate_asym_infection,
      sym_infection = rate_sym_infection,
      sym_recovery = rate_sym_recovery,
      asym_recovery = rate_asym_recovery,
      sym_removal = rate_sym_removal,
      efs_rates
    )

    total_rate <- sum(rates)
    if (!is.finite(total_rate) || total_rate <= 0) {
      break
    }

    time <- time + stats::rexp(1L, rate = total_rate)
    event <- sample.int(length(rates), size = 1L, prob = rates)

    if (event == 1L) {
      S <- S - 1L
      A <- A + 1L
    } else if (event == 2L) {
      S <- S - 1L
      I <- I + 1L
    } else if (event == 3L) {
      I <- I - 1L
      S <- S + 1L
    } else if (event == 4L) {
      A <- A - 1L
      S <- S + 1L
    } else if (event == 5L) {
      I <- I - 1L
      R <- R + 1L
    } else {
      F <- efs_targets[event - 5L]
    }

    maximum_infectious <- max(maximum_infectious, I + A)
  }

  c(
    maximum_infectious = maximum_infectious,
    extinct = as.integer(I + A == 0L),
    total_time = time,
    events = event_count
  )
}

run_gillespie_validation <- function(
    B,
    parameters,
    init,
    seed = 123L,
    max_events = 1000000L
) {
  B <- .assert_integerish(B, "B", lower = 1L)
  seed <- .assert_integerish(seed, "seed", lower = 0L)
  set.seed(seed)

  maxima <- integer(B)
  extinct <- integer(B)
  total_time <- numeric(B)
  events <- integer(B)

  for (b in seq_len(B)) {
    simulation <- simulate_variant1_gillespie(
      parameters = parameters,
      init = init,
      max_events = max_events
    )

    maxima[b] <- as.integer(simulation["maximum_infectious"])
    extinct[b] <- as.integer(simulation["extinct"])
    total_time[b] <- simulation["total_time"]
    events[b] <- as.integer(simulation["events"])
  }

  support <- 0:parameters$N
  counts <- tabulate(maxima + 1L, nbins = parameters$N + 1L)
  pmf <- counts / B
  cdf <- cumsum(pmf)

  list(
    distribution = data.frame(
      maximum = support,
      probability = pmf,
      cumulative_probability = cdf
    ),
    summary = data.frame(
      simulations = B,
      mean = mean(maxima),
      variance = stats::var(maxima),
      sd = stats::sd(maxima),
      extinction_probability = mean(extinct),
      mean_extinction_time = mean(total_time),
      mean_events = mean(events)
    ),
    simulations = data.frame(
      maximum_infectious = maxima,
      extinct = extinct,
      total_time = total_time,
      events = events
    )
  )
}

plot_maximum_comparison <- function(comparison, file) {
  grDevices::png(filename = file, width = 1600, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)

  probabilities <- rbind(
    Exact = comparison$exact_probability,
    Gillespie = comparison$gillespie_probability
  )

  graphics::barplot(
    probabilities,
    beside = TRUE,
    names.arg = comparison$maximum,
    xlab = "Maximum number of simultaneously infectious individuals",
    ylab = "Probability",
    main = "SuASR Variant I: exact and Gillespie distributions",
    legend.text = rownames(probabilities),
    args.legend = list(x = "topright", bty = "n"),
    las = 1,
    cex.names = 0.8
  )
}

run_single_example <- function(
    B = 5000L,
    seed = 123L,
    output_dir = "outputs",
    make_plot = TRUE
) {
  # Manuscript validation case 1.
  parameters <- new_variant1_parameters(
    N = 20,
    p = 0.4,
    beta = c(1.0, 1.5, 0.5),
    gamma = c(0.2, 0.2, 0.2),
    kappa = c(0.3, 0.3, 0.3),
    mu = c(0.15, 0.05, 0.01),
    Q = matrix(
      c(
        -1 / 5,  1 / 10,  1 / 10,
         1 / 15, -2 / 15,  1 / 15,
         1 / 25,  1 / 25, -2 / 25
      ),
      nrow = 3,
      byrow = TRUE
    )
  )

  initial_state <- c(S = 13, I = 1, A = 1, R = 5, F = 1)

  message("Computing the exact maximum distribution...")
  exact <- exact_maximum_distribution(
    parameters = parameters,
    init = initial_state
  )

  message(sprintf("Running %d Gillespie simulations...", B))
  gillespie <- run_gillespie_validation(
    B = B,
    parameters = parameters,
    init = initial_state,
    seed = seed
  )

  comparison <- data.frame(
    maximum = exact$distribution$maximum,
    exact_probability = exact$distribution$probability,
    gillespie_probability = gillespie$distribution$probability,
    exact_cdf = exact$distribution$cumulative_probability,
    gillespie_cdf = gillespie$distribution$cumulative_probability
  )
  comparison$absolute_probability_difference <- abs(
    comparison$exact_probability - comparison$gillespie_probability
  )
  comparison$absolute_cdf_difference <- abs(
    comparison$exact_cdf - comparison$gillespie_cdf
  )

  summary_table <- data.frame(
    method = c("Exact", "Gillespie"),
    mean = c(exact$summary$mean, gillespie$summary$mean),
    variance = c(exact$summary$variance, gillespie$summary$variance),
    sd = c(exact$summary$sd, gillespie$summary$sd),
    probability_mass_or_extinction = c(
      exact$summary$probability_mass,
      gillespie$summary$extinction_probability
    )
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(
    exact$distribution,
    file = file.path(output_dir, "variant1_maximum_exact.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    gillespie$distribution,
    file = file.path(output_dir, "variant1_maximum_gillespie.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    comparison,
    file = file.path(output_dir, "variant1_maximum_comparison.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary_table,
    file = file.path(output_dir, "variant1_maximum_summary.csv"),
    row.names = FALSE
  )

  if (isTRUE(make_plot)) {
    plot_maximum_comparison(
      comparison,
      file = file.path(output_dir, "variant1_maximum_comparison.png")
    )
  }

  cat("\nExample: Manuscript validation case 1\n")
  cat("Initial state: (S, I, A, R, F) = (13, 1, 1, 5, 1)\n\n")
  print(summary_table, row.names = FALSE)

  cat(sprintf(
    "\nMaximum absolute CDF difference: %.6f\n",
    max(comparison$absolute_cdf_difference)
  ))
  cat(sprintf("Results written to: %s\n", normalizePath(output_dir)))

  invisible(
    list(
      parameters = parameters,
      initial_state = initial_state,
      exact = exact,
      gillespie = gillespie,
      comparison = comparison,
      summary = summary_table
    )
  )
}

if (sys.nframe() == 0L) {
  run_single_example()
}
