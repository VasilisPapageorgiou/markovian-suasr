# SuASR Variant I: exact joint infection-count distribution
#
# This file implements the recurrent-infection formulation of the Markovian
# SuASR model with stochastic external factor states. It computes the joint
# distribution of the numbers of additional symptomatic and asymptomatic
# infection events until extinction by the matrix recursion in Equation (19)
# of the associated manuscript.
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

# Standardize an external factor generator.
# Only the off-diagonal entries supplied by the user are used. The diagonal is
# recalculated as minus the sum of the off-diagonal entries in each row.
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

# Create and validate a parameter object for Variant I.
new_suasr_variant1_parameters <- function(
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
  invalid_rate <- vapply(rate_vectors, function(x) any(x < 0), logical(1))
  if (any(invalid_rate)) {
    stop(
      sprintf(
        "The following rate vectors contain negative values: %s.",
        paste(names(rate_vectors)[invalid_rate], collapse = ", ")
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

print.suasr_variant1_parameters <- function(x, ...) {
  cat("SuASR Variant I parameters\n")
  cat(sprintf("  Population size: %d\n", x$N))
  cat(sprintf("  External factor states: %d\n", x$n_efs))
  invisible(x)
}

.validate_initial_state <- function(init, parameters, require_infectious = TRUE) {
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

  if (require_infectious && I + A == 0L) {
    stop("The initial state must contain at least one infectious individual.", call. = FALSE)
  }

  c(S = S, I = I, A = A, R = R, F = F)
}

.enumerate_transient_states <- function(parameters) {
  N <- parameters$N
  n_base <- as.integer(round(choose(N + 3L, 3L)))
  base <- matrix(0L, nrow = n_base, ncol = 4L)
  colnames(base) <- c("S", "I", "A", "R")

  position <- 1L
  for (I in 0:N) {
    for (A in 0:(N - I)) {
      for (R in 0:(N - I - A)) {
        S <- N - I - A - R
        base[position, ] <- c(S, I, A, R)
        position <- position + 1L
      }
    }
  }

  base <- as.data.frame(base, stringsAsFactors = FALSE)
  base <- base[base$I + base$A > 0L, , drop = FALSE]

  states <- do.call(
    rbind,
    lapply(seq_len(parameters$n_efs), function(F) {
      cbind(base, F = F)
    })
  )
  rownames(states) <- NULL

  integer_columns <- c("S", "I", "A", "R", "F")
  states[integer_columns] <- lapply(states[integer_columns], as.integer)
  states
}

.state_key <- function(I, A, R, F) {
  paste(I, A, R, F, sep = ":")
}

# Build the sparse matrices used in Equation (19).
#
# U_core: transitions that do not create a new infection event, with a
#         diagonal containing the negative total rate of every possible event.
# U_sym:  transitions that create one new symptomatic infection.
# U_asym: transitions that create one new asymptomatic infection.
# q_abs:  rates from transient states to the disease-free set.
build_variant1_generator_blocks <- function(parameters, check = TRUE) {
  .require_matrix()

  if (!inherits(parameters, "suasr_variant1_parameters")) {
    stop(
      "'parameters' must be created by new_suasr_variant1_parameters().",
      call. = FALSE
    )
  }

  states <- .enumerate_transient_states(parameters)
  n_states <- nrow(states)

  keys <- .state_key(states$I, states$A, states$R, states$F)
  state_lookup <- setNames(seq_len(n_states), keys)

  lookup_state <- function(I, A, R, F) {
    key <- .state_key(I, A, R, F)
    id <- unname(state_lookup[key])
    if (length(id) != 1L || is.na(id)) {
      stop(sprintf("Internal error: state %s was not found.", key), call. = FALSE)
    }
    as.integer(id)
  }

  core_row_list <- vector("list", n_states)
  core_col_list <- vector("list", n_states)
  core_value_list <- vector("list", n_states)

  sym_rows <- integer(n_states)
  sym_cols <- integer(n_states)
  sym_values <- numeric(n_states)
  sym_count <- 0L

  asym_rows <- integer(n_states)
  asym_cols <- integer(n_states)
  asym_values <- numeric(n_states)
  asym_count <- 0L

  diagonal <- numeric(n_states)
  q_abs <- numeric(n_states)

  for (row in seq_len(n_states)) {
    S <- states$S[row]
    I <- states$I[row]
    A <- states$A[row]
    R <- states$R[row]
    F <- states$F[row]

    row_cols <- integer(0)
    row_values <- numeric(0)
    total_rate <- 0

    infection_force <- I + A
    base_infection_rate <- parameters$beta[F] * S * infection_force / parameters$N
    symptomatic_rate <- parameters$p[F] * base_infection_rate
    asymptomatic_rate <- (1 - parameters$p[F]) * base_infection_rate

    total_rate <- total_rate + symptomatic_rate + asymptomatic_rate

    if (symptomatic_rate > 0) {
      sym_count <- sym_count + 1L
      sym_rows[sym_count] <- row
      sym_cols[sym_count] <- lookup_state(I + 1L, A, R, F)
      sym_values[sym_count] <- symptomatic_rate
    }

    if (asymptomatic_rate > 0) {
      asym_count <- asym_count + 1L
      asym_rows[asym_count] <- row
      asym_cols[asym_count] <- lookup_state(I, A + 1L, R, F)
      asym_values[asym_count] <- asymptomatic_rate
    }

    if (I > 0L) {
      symptomatic_recovery_rate <- parameters$gamma[F] * I
      total_rate <- total_rate + symptomatic_recovery_rate

      if (symptomatic_recovery_rate > 0) {
        if ((I - 1L) + A == 0L) {
          q_abs[row] <- q_abs[row] + symptomatic_recovery_rate
        } else {
          row_cols <- c(row_cols, lookup_state(I - 1L, A, R, F))
          row_values <- c(row_values, symptomatic_recovery_rate)
        }
      }

      symptomatic_removal_rate <- parameters$mu[F] * I
      total_rate <- total_rate + symptomatic_removal_rate

      if (symptomatic_removal_rate > 0) {
        if ((I - 1L) + A == 0L) {
          q_abs[row] <- q_abs[row] + symptomatic_removal_rate
        } else {
          row_cols <- c(row_cols, lookup_state(I - 1L, A, R + 1L, F))
          row_values <- c(row_values, symptomatic_removal_rate)
        }
      }
    }

    if (A > 0L) {
      asymptomatic_recovery_rate <- parameters$kappa[F] * A
      total_rate <- total_rate + asymptomatic_recovery_rate

      if (asymptomatic_recovery_rate > 0) {
        if (I + (A - 1L) == 0L) {
          q_abs[row] <- q_abs[row] + asymptomatic_recovery_rate
        } else {
          row_cols <- c(row_cols, lookup_state(I, A - 1L, R, F))
          row_values <- c(row_values, asymptomatic_recovery_rate)
        }
      }
    }

    efs_targets <- which(parameters$Q[F, ] > 0)
    if (length(efs_targets) > 0L) {
      for (target_F in efs_targets) {
        efs_rate <- parameters$Q[F, target_F]
        total_rate <- total_rate + efs_rate
        row_cols <- c(row_cols, lookup_state(I, A, R, target_F))
        row_values <- c(row_values, efs_rate)
      }
    }

    diagonal[row] <- -total_rate
    core_row_list[[row]] <- rep.int(row, length(row_cols))
    core_col_list[[row]] <- row_cols
    core_value_list[[row]] <- row_values
  }

  core_rows <- c(seq_len(n_states), unlist(core_row_list, use.names = FALSE))
  core_cols <- c(seq_len(n_states), unlist(core_col_list, use.names = FALSE))
  core_values <- c(diagonal, unlist(core_value_list, use.names = FALSE))

  U_core <- Matrix::sparseMatrix(
    i = core_rows,
    j = core_cols,
    x = core_values,
    dims = c(n_states, n_states),
    giveCsparse = TRUE
  )

  if (sym_count > 0L) {
    U_sym <- Matrix::sparseMatrix(
      i = sym_rows[seq_len(sym_count)],
      j = sym_cols[seq_len(sym_count)],
      x = sym_values[seq_len(sym_count)],
      dims = c(n_states, n_states),
      giveCsparse = TRUE
    )
  } else {
    U_sym <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_states, n_states), giveCsparse = TRUE
    )
  }

  if (asym_count > 0L) {
    U_asym <- Matrix::sparseMatrix(
      i = asym_rows[seq_len(asym_count)],
      j = asym_cols[seq_len(asym_count)],
      x = asym_values[seq_len(asym_count)],
      dims = c(n_states, n_states),
      giveCsparse = TRUE
    )
  } else {
    U_asym <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_states, n_states), giveCsparse = TRUE
    )
  }

  if (check) {
    balance <- as.numeric(rowSums(U_core)) +
      as.numeric(rowSums(U_sym)) +
      as.numeric(rowSums(U_asym)) + q_abs

    maximum_error <- max(abs(balance))
    if (!is.finite(maximum_error) || maximum_error > 1e-9) {
      stop(
        sprintf(
          "Generator balance check failed; maximum absolute row error is %.3e.",
          maximum_error
        ),
        call. = FALSE
      )
    }
  }

  structure(
    list(
      parameters = parameters,
      states = states,
      state_lookup = state_lookup,
      U_core = U_core,
      U_sym = U_sym,
      U_asym = U_asym,
      q_abs = q_abs
    ),
    class = "suasr_variant1_blocks"
  )
}

print.suasr_variant1_blocks <- function(x, ...) {
  cat("SuASR Variant I generator blocks\n")
  cat(sprintf("  Transient states: %d\n", nrow(x$states)))
  cat(sprintf("  Nonzero entries in U_core: %d\n", length(x$U_core@x)))
  cat(sprintf("  Nonzero entries in U_sym: %d\n", length(x$U_sym@x)))
  cat(sprintf("  Nonzero entries in U_asym: %d\n", length(x$U_asym@x)))
  invisible(x)
}

# Compute the truncated joint distribution for one selected initial state.
#
# Rows correspond to 0:max_sym additional symptomatic infection events.
# Columns correspond to 0:max_asym additional asymptomatic infection events.
#
# The anti-diagonal implementation stores only the previous and current
# anti-diagonals of probability vectors, rather than the full three-dimensional
# tensor. This substantially reduces memory use while preserving Equation (19).
compute_variant1_joint_distribution <- function(
    blocks,
    init,
    max_sym = 100L,
    max_asym = 100L,
    progress = interactive(),
    negative_tolerance = 1e-10
) {
  .require_matrix()

  if (!inherits(blocks, "suasr_variant1_blocks")) {
    stop(
      "'blocks' must be created by build_variant1_generator_blocks().",
      call. = FALSE
    )
  }

  max_sym <- .assert_integerish(max_sym, "max_sym", lower = 0L)
  max_asym <- .assert_integerish(max_asym, "max_asym", lower = 0L)
  init <- .validate_initial_state(init, blocks$parameters, require_infectious = TRUE)

  initial_key <- .state_key(init["I"], init["A"], init["R"], init["F"])
  initial_id <- unname(blocks$state_lookup[initial_key])
  if (length(initial_id) != 1L || is.na(initial_id)) {
    stop("The selected initial state is not in the transient state space.", call. = FALSE)
  }
  initial_id <- as.integer(initial_id)

  lu_factor <- tryCatch(
    Matrix::lu(blocks$U_core),
    error = function(error) {
      stop(
        paste0(
          "Sparse LU factorization failed. Check whether the transient block ",
          "is nonsingular and absorption occurs with probability one. Original error: ",
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  solve_core <- function(rhs) {
    solution <- tryCatch(
      solve(lu_factor, rhs),
      error = function(error) {
        stop(
          paste0("Linear solve failed: ", conditionMessage(error)),
          call. = FALSE
        )
      }
    )
    as.matrix(solution)
  }

  n_states <- nrow(blocks$states)
  joint <- matrix(
    0,
    nrow = max_sym + 1L,
    ncol = max_asym + 1L,
    dimnames = list(
      N_symptomatic = as.character(0:max_sym),
      N_asymptomatic = as.character(0:max_asym)
    )
  )

  p00 <- -as.numeric(solve_core(matrix(blocks$q_abs, ncol = 1L)))
  previous_values <- matrix(p00, nrow = n_states, ncol = 1L)
  previous_n <- 0L
  joint[1L, 1L] <- p00[initial_id]

  maximum_diagonal <- max_sym + max_asym
  progress_bar <- NULL
  if (isTRUE(progress) && maximum_diagonal > 0L) {
    progress_bar <- utils::txtProgressBar(
      min = 0L,
      max = maximum_diagonal,
      style = 3L
    )
    on.exit(close(progress_bar), add = TRUE)
  }

  if (maximum_diagonal > 0L) {
    for (diagonal_index in seq_len(maximum_diagonal)) {
      n_values <- seq.int(
        from = max(0L, diagonal_index - max_asym),
        to = min(max_sym, diagonal_index)
      )
      m_values <- diagonal_index - n_values
      n_columns <- length(n_values)

      rhs <- matrix(0, nrow = n_states, ncol = n_columns)

      has_symptomatic_predecessor <- n_values > 0L
      if (any(has_symptomatic_predecessor)) {
        predecessor_columns <- match(
          n_values[has_symptomatic_predecessor] - 1L,
          previous_n
        )
        if (anyNA(predecessor_columns)) {
          stop("Internal error in symptomatic anti-diagonal indexing.", call. = FALSE)
        }

        rhs[, has_symptomatic_predecessor] <- as.matrix(
          blocks$U_sym %*%
            previous_values[, predecessor_columns, drop = FALSE]
        )
      }

      has_asymptomatic_predecessor <- m_values > 0L
      if (any(has_asymptomatic_predecessor)) {
        predecessor_columns <- match(
          n_values[has_asymptomatic_predecessor],
          previous_n
        )
        if (anyNA(predecessor_columns)) {
          stop("Internal error in asymptomatic anti-diagonal indexing.", call. = FALSE)
        }

        rhs[, has_asymptomatic_predecessor] <-
          rhs[, has_asymptomatic_predecessor, drop = FALSE] +
          as.matrix(
            blocks$U_asym %*%
              previous_values[, predecessor_columns, drop = FALSE]
          )
      }

      current_values <- -solve_core(rhs)
      if (ncol(current_values) != n_columns) {
        current_values <- matrix(
          as.numeric(current_values),
          nrow = n_states,
          ncol = n_columns
        )
      }

      joint[cbind(n_values + 1L, m_values + 1L)] <-
        current_values[initial_id, ]

      previous_values <- current_values
      previous_n <- n_values

      if (!is.null(progress_bar)) {
        utils::setTxtProgressBar(progress_bar, diagonal_index)
      }
    }
  }

  minimum_probability <- min(joint)
  if (minimum_probability < -negative_tolerance) {
    warning(
      sprintf(
        "The computed grid contains a negative probability as small as %.3e.",
        minimum_probability
      ),
      call. = FALSE
    )
  }
  joint[joint < 0] <- 0

  captured_mass <- sum(joint)
  omitted_mass <- 1 - captured_mass

  if (captured_mass > 1 + 1e-8) {
    warning(
      sprintf("Captured probability mass exceeds one: %.12f.", captured_mass),
      call. = FALSE
    )
  }

  symptomatic_values <- 0:max_sym
  asymptomatic_values <- 0:max_asym
  symptomatic_marginal <- rowSums(joint)
  asymptomatic_marginal <- colSums(joint)

  if (captured_mass > 0) {
    mean_symptomatic <- sum(symptomatic_values * symptomatic_marginal) / captured_mass
    mean_asymptomatic <- sum(asymptomatic_values * asymptomatic_marginal) / captured_mass

    sd_symptomatic <- sqrt(
      sum((symptomatic_values - mean_symptomatic)^2 * symptomatic_marginal) /
        captured_mass
    )
    sd_asymptomatic <- sqrt(
      sum((asymptomatic_values - mean_asymptomatic)^2 * asymptomatic_marginal) /
        captured_mass
    )
  } else {
    mean_symptomatic <- NA_real_
    mean_asymptomatic <- NA_real_
    sd_symptomatic <- NA_real_
    sd_asymptomatic <- NA_real_
  }

  summary <- data.frame(
    captured_mass = captured_mass,
    omitted_mass = omitted_mass,
    mean_symptomatic_conditional_on_grid = mean_symptomatic,
    sd_symptomatic_conditional_on_grid = sd_symptomatic,
    mean_asymptomatic_conditional_on_grid = mean_asymptomatic,
    sd_asymptomatic_conditional_on_grid = sd_asymptomatic,
    row.names = NULL
  )

  structure(
    list(
      probabilities = joint,
      symptomatic_marginal = symptomatic_marginal,
      asymptomatic_marginal = asymptomatic_marginal,
      summary = summary,
      initial_state = init,
      max_sym = max_sym,
      max_asym = max_asym,
      captured_first_moment_symptomatic =
        sum(symptomatic_values * symptomatic_marginal),
      captured_first_moment_asymptomatic =
        sum(asymptomatic_values * asymptomatic_marginal)
    ),
    class = "suasr_variant1_joint_distribution"
  )
}

print.suasr_variant1_joint_distribution <- function(x, ...) {
  cat("SuASR Variant I joint infection-count distribution\n")
  cat(
    sprintf(
      "  Initial state: S=%d, I=%d, A=%d, R=%d, F=%d\n",
      x$initial_state["S"],
      x$initial_state["I"],
      x$initial_state["A"],
      x$initial_state["R"],
      x$initial_state["F"]
    )
  )
  cat(sprintf("  Grid: symptomatic 0:%d, asymptomatic 0:%d\n", x$max_sym, x$max_asym))
  print(x$summary, row.names = FALSE, digits = 7)
  invisible(x)
}

joint_distribution_to_long <- function(result, drop_zero = TRUE) {
  if (!inherits(result, "suasr_variant1_joint_distribution")) {
    stop(
      "'result' must be returned by compute_variant1_joint_distribution().",
      call. = FALSE
    )
  }

  output <- expand.grid(
    N_symptomatic = 0:result$max_sym,
    N_asymptomatic = 0:result$max_asym,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  output$probability <- as.vector(result$probabilities)

  if (isTRUE(drop_zero)) {
    output <- output[output$probability > 0, , drop = FALSE]
    rownames(output) <- NULL
  }

  output
}

write_joint_distribution_csv <- function(result, file, drop_zero = FALSE) {
  output <- joint_distribution_to_long(result, drop_zero = drop_zero)
  utils::write.csv(output, file = file, row.names = FALSE)
  invisible(normalizePath(file, mustWork = FALSE))
}

plot_joint_distribution <- function(result) {
  if (!inherits(result, "suasr_variant1_joint_distribution")) {
    stop(
      "'result' must be returned by compute_variant1_joint_distribution().",
      call. = FALSE
    )
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "Package 'plotly' is required for this plot. Install it with install.packages('plotly').",
      call. = FALSE
    )
  }

  figure <- plotly::plot_ly(
    x = 0:result$max_sym,
    y = 0:result$max_asym,
    z = t(result$probabilities),
    type = "surface"
  )

  plotly::layout(
    figure,
    scene = list(
      xaxis = list(title = "Additional symptomatic infections"),
      yaxis = list(title = "Additional asymptomatic infections"),
      zaxis = list(title = "Probability")
    )
  )
}

# One exact Gillespie trajectory for Variant I.
simulate_variant1_gillespie_once <- function(
    parameters,
    init,
    max_events = 1000000L
) {
  if (!inherits(parameters, "suasr_variant1_parameters")) {
    stop(
      "'parameters' must be created by new_suasr_variant1_parameters().",
      call. = FALSE
    )
  }

  max_events <- .assert_integerish(max_events, "max_events", lower = 1L)
  state <- .validate_initial_state(init, parameters, require_infectious = FALSE)

  S <- state["S"]
  I <- state["I"]
  A <- state["A"]
  R <- state["R"]
  F <- state["F"]

  time <- 0
  symptomatic_count <- 0L
  asymptomatic_count <- 0L
  maximum_infectious <- I + A
  event_count <- 0L

  while (I + A > 0L && event_count < max_events) {
    infection_force <- I + A
    base_infection_rate <- parameters$beta[F] * S * infection_force / parameters$N

    fixed_rates <- c(
      asymptomatic_infection = (1 - parameters$p[F]) * base_infection_rate,
      symptomatic_infection = parameters$p[F] * base_infection_rate,
      symptomatic_recovery = parameters$gamma[F] * I,
      asymptomatic_recovery = parameters$kappa[F] * A,
      symptomatic_removal = parameters$mu[F] * I
    )

    efs_targets <- which(parameters$Q[F, ] > 0)
    efs_rates <- parameters$Q[F, efs_targets]
    rates <- c(fixed_rates, efs_rates)

    total_rate <- sum(rates)
    if (!is.finite(total_rate) || total_rate <= 0) {
      break
    }

    time <- time + stats::rexp(1L, rate = total_rate)
    selected_event <- sample.int(length(rates), size = 1L, prob = rates)
    event_count <- event_count + 1L

    if (selected_event == 1L) {
      S <- S - 1L
      A <- A + 1L
      asymptomatic_count <- asymptomatic_count + 1L
    } else if (selected_event == 2L) {
      S <- S - 1L
      I <- I + 1L
      symptomatic_count <- symptomatic_count + 1L
    } else if (selected_event == 3L) {
      I <- I - 1L
      S <- S + 1L
    } else if (selected_event == 4L) {
      A <- A - 1L
      S <- S + 1L
    } else if (selected_event == 5L) {
      I <- I - 1L
      R <- R + 1L
    } else {
      F <- efs_targets[selected_event - 5L]
    }

    maximum_infectious <- max(maximum_infectious, I + A)
  }

  data.frame(
    total_time = time,
    N_symptomatic_additional = symptomatic_count,
    N_asymptomatic_additional = asymptomatic_count,
    N_total_additional = symptomatic_count + asymptomatic_count,
    maximum_infectious = maximum_infectious,
    extinct = as.integer(I + A == 0L),
    stopped_at_max_events = as.integer(I + A > 0L && event_count >= max_events),
    final_S = S,
    final_I = I,
    final_A = A,
    final_R = R,
    final_F = F,
    events = event_count,
    row.names = NULL
  )
}

run_variant1_gillespie <- function(
    parameters,
    init,
    n_simulations = 5000L,
    seed = 123L,
    max_events = 1000000L,
    progress = interactive()
) {
  n_simulations <- .assert_integerish(
    n_simulations,
    "n_simulations",
    lower = 1L
  )
  max_events <- .assert_integerish(max_events, "max_events", lower = 1L)

  if (!is.null(seed)) {
    seed <- .assert_integerish(seed, "seed", lower = 0L)
    set.seed(seed)
  }

  simulations <- vector("list", n_simulations)
  progress_bar <- NULL
  if (isTRUE(progress)) {
    progress_bar <- utils::txtProgressBar(
      min = 0L,
      max = n_simulations,
      style = 3L
    )
    on.exit(close(progress_bar), add = TRUE)
  }

  for (simulation_index in seq_len(n_simulations)) {
    simulations[[simulation_index]] <- simulate_variant1_gillespie_once(
      parameters = parameters,
      init = init,
      max_events = max_events
    )

    if (!is.null(progress_bar)) {
      utils::setTxtProgressBar(progress_bar, simulation_index)
    }
  }

  simulations <- do.call(rbind, simulations)

  if (any(simulations$stopped_at_max_events == 1L)) {
    warning(
      "At least one trajectory reached 'max_events' before extinction.",
      call. = FALSE
    )
  }

  summary <- data.frame(
    n_simulations = n_simulations,
    mean_symptomatic = mean(simulations$N_symptomatic_additional),
    sd_symptomatic = stats::sd(simulations$N_symptomatic_additional),
    mean_asymptomatic = mean(simulations$N_asymptomatic_additional),
    sd_asymptomatic = stats::sd(simulations$N_asymptomatic_additional),
    mean_total = mean(simulations$N_total_additional),
    sd_total = stats::sd(simulations$N_total_additional),
    mean_maximum_infectious = mean(simulations$maximum_infectious),
    sd_maximum_infectious = stats::sd(simulations$maximum_infectious),
    extinction_fraction = mean(simulations$extinct),
    stopped_at_max_events_fraction = mean(simulations$stopped_at_max_events),
    row.names = NULL
  )

  joint_table <- table(
    N_symptomatic = simulations$N_symptomatic_additional,
    N_asymptomatic = simulations$N_asymptomatic_additional
  )
  joint_distribution <- as.data.frame(
    prop.table(joint_table),
    responseName = "probability",
    stringsAsFactors = FALSE
  )
  joint_distribution$N_symptomatic <-
    as.integer(as.character(joint_distribution$N_symptomatic))
  joint_distribution$N_asymptomatic <-
    as.integer(as.character(joint_distribution$N_asymptomatic))
  joint_distribution <- joint_distribution[
    joint_distribution$probability > 0,
    ,
    drop = FALSE
  ]
  rownames(joint_distribution) <- NULL

  maximum_table <- table(simulations$maximum_infectious)
  maximum_distribution <- data.frame(
    maximum_infectious = as.integer(names(maximum_table)),
    probability = as.numeric(prop.table(maximum_table)),
    row.names = NULL
  )

  structure(
    list(
      simulations = simulations,
      summary = summary,
      joint_distribution = joint_distribution,
      maximum_distribution = maximum_distribution
    ),
    class = "suasr_variant1_gillespie"
  )
}

print.suasr_variant1_gillespie <- function(x, ...) {
  cat("SuASR Variant I Gillespie results\n")
  print(x$summary, row.names = FALSE, digits = 7)
  invisible(x)
}

# ============================================================================
# Single reproducible example: manuscript validation case 1
# ============================================================================
# Run from the command line with:
#   Rscript suasr_variant1_single_example.R
#
# When the file is sourced, call:
#   run_single_example()

run_single_example <- function(
    output_directory = "results",
    max_symptomatic = 149L,
    max_asymptomatic = 149L,
    n_gillespie = 5000L,
    seed = 123L,
    progress = interactive()
) {
  dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

  # External-factor generator. Diagonal entries are recalculated internally.
  Q <- matrix(
    c(
      0,      1 / 10, 1 / 10,
      1 / 15, 0,      1 / 15,
      1 / 25, 1 / 25, 0
    ),
    nrow = 3L,
    byrow = TRUE
  )

  parameters <- new_suasr_variant1_parameters(
    N = 20L,
    p = 0.4,
    beta = c(1.0, 1.5, 0.5),
    gamma = c(0.2, 0.2, 0.2),
    kappa = c(0.3, 0.3, 0.3),
    mu = c(0.15, 0.05, 0.01),
    Q = Q
  )

  # Counts generated after time zero exclude the initially infectious subjects.
  initial_state <- c(S = 13L, I = 1L, A = 1L, R = 5L, F = 1L)

  message("Building sparse generator blocks...")
  blocks <- build_variant1_generator_blocks(parameters)
  print(blocks)

  message("Computing the exact truncated joint distribution...")
  exact_result <- compute_variant1_joint_distribution(
    blocks = blocks,
    init = initial_state,
    max_sym = max_symptomatic,
    max_asym = max_asymptomatic,
    progress = progress
  )
  print(exact_result)

  exact_distribution_file <- file.path(
    output_directory,
    "exact_joint_distribution.csv"
  )
  exact_summary_file <- file.path(
    output_directory,
    "exact_summary.csv"
  )

  write_joint_distribution_csv(
    result = exact_result,
    file = exact_distribution_file,
    drop_zero = FALSE
  )
  utils::write.csv(
    exact_result$summary,
    exact_summary_file,
    row.names = FALSE
  )

  message("Running Gillespie validation simulations...")
  gillespie_result <- run_variant1_gillespie(
    parameters = parameters,
    init = initial_state,
    n_simulations = n_gillespie,
    seed = seed,
    progress = progress
  )
  print(gillespie_result)

  utils::write.csv(
    gillespie_result$simulations,
    file.path(output_directory, "gillespie_trajectories.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    gillespie_result$joint_distribution,
    file.path(output_directory, "gillespie_joint_distribution.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    gillespie_result$summary,
    file.path(output_directory, "gillespie_summary.csv"),
    row.names = FALSE
  )

  exact_summary <- exact_result$summary
  simulation_summary <- gillespie_result$summary

  comparison <- data.frame(
    method = c("Exact matrix recursion", "Gillespie simulation"),
    mean_symptomatic = c(
      exact_summary$mean_symptomatic_conditional_on_grid,
      simulation_summary$mean_symptomatic
    ),
    sd_symptomatic = c(
      exact_summary$sd_symptomatic_conditional_on_grid,
      simulation_summary$sd_symptomatic
    ),
    mean_asymptomatic = c(
      exact_summary$mean_asymptomatic_conditional_on_grid,
      simulation_summary$mean_asymptomatic
    ),
    sd_asymptomatic = c(
      exact_summary$sd_asymptomatic_conditional_on_grid,
      simulation_summary$sd_asymptomatic
    ),
    captured_mass = c(exact_summary$captured_mass, NA_real_),
    row.names = NULL
  )

  utils::write.csv(
    comparison,
    file.path(output_directory, "exact_vs_gillespie_summary.csv"),
    row.names = FALSE
  )

  message("Completed. Output files were written to: ", output_directory)
  print(comparison, row.names = FALSE)

  invisible(
    list(
      parameters = parameters,
      initial_state = initial_state,
      blocks = blocks,
      exact = exact_result,
      gillespie = gillespie_result,
      comparison = comparison
    )
  )
}

if (sys.nframe() == 0L) {
  run_single_example()
}
