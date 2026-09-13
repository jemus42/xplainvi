#' @title Shapley Additive Global Importance (SAGE) Base Class
#'
#' @description Base class for SAGE (Shapley Additive Global Importance)
#' feature importance based on Shapley values with marginalization.
#' This is an abstract class - use [MarginalSAGE] or [ConditionalSAGE].
#'
#' @details
#' SAGE uses Shapley values to fairly distribute the total prediction
#' performance among all features. Unlike perturbation-based methods,
#' SAGE marginalizes features by integrating over their distribution.
#' This is approximated by averaging predictions over a reference dataset.
#'
#' SAGE values are reductions in the measure's score relative to the empty coalition,
#' `score(empty) - score(S)`, so that positive values mean the feature improves performance.
#' For measures that are maximized (`measure$minimize = FALSE`, e.g. `classif.acc`) the scores are
#' negated internally, so the sign convention is the same for all measures.
#'
#' **Standard Error Calculation**: The standard errors (SE) reported in
#' `$convergence_history` reflect the uncertainty in Shapley value estimation
#' across different random permutations within a single resampling iteration.
#' These SEs quantify the Monte Carlo sampling error for a fixed trained model
#' and are only valid for inference about the importance of features for that
#' specific model. They do not capture broader uncertainty from model variability
#' across different train/test splits or resampling iterations.
#'
#' @references
#' `r print_bib("lundberg_2020")`
#'
#' @seealso [MarginalSAGE] [ConditionalSAGE]
#'
#' @export
SAGE = R6Class(
  "SAGE",
  inherit = FeatureImportanceMethod,
  public = list(
    #' @field n_permutations (`integer(1)`) Number of permutations to sample.
    n_permutations = NULL,
    #' @field convergence_history ([`data.table`][data.table::data.table]) History of SAGE values during computation.
    convergence_history = NULL,
    #' @field converged (`logical(1)`) Whether convergence was detected.
    converged = FALSE,
    #' @field n_permutations_used (`integer(1)`) Actual number of permutations used.
    n_permutations_used = NULL,

    #' @description
    #' Creates a new instance of the SAGE class.
    #' @param task,learner,measure,resampling,features Passed to FeatureImportanceMethod.
    #' @param n_permutations (`integer(1)`: `10L`) Number of permutations to sample for SAGE value estimation.
    #'   The total number of evaluated coalitions is `1 (empty) + n_permutations * n_features`.
    #' @param batch_size (`integer(1)`: `5000L`) Maximum number of observations to process in a single prediction call.
    #' @param n_samples (`integer(1)`: `100L`) Number of samples to use for marginalizing out-of-coalition features.
    #'   For [MarginalSAGE], this is the number of marginal data samples ("background data" in other implementations).
    #'   For [ConditionalSAGE], this is the number of conditional samples per test instance retrieved from `sampler`.
    #' @param early_stopping (`logical(1)`: `TRUE`) Whether to enable early stopping based on convergence detection.
    #' @param se_threshold (`numeric(1)`: `0.01`) Convergence threshold for relative standard error.
    #'   Convergence is detected when the maximum relative SE across all features falls below this threshold.
    #'   Relative SE is calculated as SE divided by the range of importance values (max - min),
    #'   making it scale-invariant across different loss metrics.
    #'   Default of `0.01` means convergence when relative SE is below 1% of the importance range.
    #' @param min_permutations (`integer(1)`: `10L`) Minimum permutations before checking for convergence. Convergence is judged based on the standard errors of the estimated SAGE values,
    #' which requires a sufficiently large number of samples (i.e., evaluated coalitions).
    #' @param check_interval (`integer(1)`: `1L`) Check convergence every N permutations.
    initialize = function(
      task,
      learner,
      measure = NULL,
      resampling = NULL,
      features = NULL,
      n_permutations = 10L,
      batch_size = 5000L,
      n_samples = 100L,
      early_stopping = TRUE,
      se_threshold = 0.01,
      min_permutations = 10L,
      check_interval = 1L
    ) {
      super$initialize(
        task = task,
        learner = learner,
        measure = measure,
        resampling = resampling,
        features = features,
        label = "Shapley Additive Global Importance"
      )

      self$n_permutations = checkmate::assert_int(n_permutations, lower = 1L)

      # For classification tasks, require predict_type = "prob"
      if (self$task$task_type == "classif") {
        if (learner$predict_type != "prob") {
          cli::cli_abort(c(
            "Classification learners require probability predictions for SAGE.",
            "i" = "Please set {.code learner$configure(predict_type = \"prob\")} before using SAGE."
          ))
        }
      }

      # Set parameters
      ps = ps(
        n_permutations = paradox::p_int(lower = 1L, default = 10L),
        batch_size = paradox::p_int(lower = 1L, default = 5000L),
        n_samples = paradox::p_int(lower = 1L, default = 100L),
        early_stopping = paradox::p_lgl(default = TRUE),
        se_threshold = paradox::p_dbl(lower = 0, upper = 1, default = 0.01),
        min_permutations = paradox::p_int(lower = 1L, default = 3L),
        check_interval = paradox::p_int(lower = 1L, default = 1L)
      )
      ps$values$n_permutations = n_permutations
      ps$values$batch_size = batch_size
      ps$values$n_samples = n_samples
      ps$values$early_stopping = early_stopping
      ps$values$se_threshold = se_threshold
      ps$values$min_permutations = min_permutations
      ps$values$check_interval = check_interval
      self$param_set = ps
    },

    #' @description
    #' Compute SAGE values.
    #' @param store_backends (`logical(1)`) Whether to store data backends.
    #' @param batch_size (`integer(1)`: `5000L`) Maximum number of observations to process in a single prediction call.
    #' @param early_stopping (`logical(1)`: `TRUE`) Whether to check for convergence and stop early.
    #' @param se_threshold (`numeric(1)`: `0.01`) Convergence threshold for relative standard error.
    #'   SE is normalized by the range of importance values (max - min) to make convergence
    #'   detection scale-invariant. Default `0.01` means convergence when relative SE < 1%.
    #' @param min_permutations (`integer(1)`: `10L`) Minimum permutations before checking convergence.
    #' @param check_interval (`integer(1)`: `1L`) Check convergence every N permutations.
    compute = function(
      store_backends = TRUE,
      batch_size = NULL,
      early_stopping = NULL,
      se_threshold = NULL,
      min_permutations = NULL,
      check_interval = NULL
    ) {
      # Reset convergence tracking
      self$convergence_history = NULL
      self$converged = FALSE
      self$n_permutations_used = NULL

      # Resolve parameters using hierarchical resolution
      batch_size = resolve_param(batch_size, self$param_set$values$batch_size, 5000L)
      early_stopping = resolve_param(
        early_stopping,
        self$param_set$values$early_stopping,
        TRUE
      )
      se_threshold = resolve_param(
        se_threshold,
        self$param_set$values$se_threshold,
        0.01
      )
      min_permutations = resolve_param(
        min_permutations,
        self$param_set$values$min_permutations,
        3L
      )
      check_interval = resolve_param(check_interval, self$param_set$values$check_interval, 1L)

      # Initial resampling to get trained learners
      rr = assemble_rr(
        task = self$task,
        learner = self$learner,
        resampling = self$resampling,
        store_models = TRUE,
        store_backends = store_backends
      )
      # Store results
      self$resample_result = rr

      # For convergence tracking, we'll use the first resampling iteration
      # (convergence is about permutation count, not resampling)
      iter_for_convergence = 1L

      # Compute SAGE values for convergence tracking (first iteration)
      first_result = private$.compute_sage_scores(
        learner = rr$learners[[iter_for_convergence]],
        test_dt = self$task$data(rows = rr$resampling$test_set(iter_for_convergence)),
        n_permutations = self$n_permutations,
        batch_size = batch_size,
        early_stopping = early_stopping,
        se_threshold = se_threshold,
        min_permutations = min_permutations,
        check_interval = check_interval
      )

      # Extract convergence data from first iteration
      # `convergence_data` exists even if early_stopping = FALSE
      self$convergence_history = first_result$convergence_data$convergence_history
      self$converged = first_result$convergence_data$converged
      self$n_permutations_used = first_result$convergence_data$n_permutations_used

      # If we have multiple resampling iterations, compute the rest without convergence tracking
      if (self$resampling$iters > 1) {
        remaining_results = lapply(seq_len(self$resampling$iters)[-iter_for_convergence], \(iter) {
          private$.compute_sage_scores(
            learner = rr$learners[[iter]],
            test_dt = self$task$data(rows = rr$resampling$test_set(iter)),
            # n_permutations_used is either same as n_permutations (if early_stopping = FALSE)
            # or its a smaller value if early_stopping = TRUE and it stopped early
            # The fallback to self$n_permutations should not be needed
            n_permutations = self$n_permutations_used %||% self$n_permutations,
            batch_size = batch_size,
            # Only track convergence etc. for first iteration
            early_stopping = FALSE
          )
        })

        # Extract scores from all results (always list format now)
        all_scores = c(list(first_result$scores), lapply(remaining_results, function(x) x$scores))
      } else {
        all_scores = list(first_result$scores)
      }

      # Combine results across resampling iterations
      scores = rbindlist(all_scores, idcol = "iter_rsmp")

      # iter_rsmp, feature, importance -- score_baseline or so don't apply here
      private$.scores = scores
    },

    #' @description
    #' Plot convergence history of SAGE values.
    #' @param features (`character` | `NULL`) Features to plot. If NULL, plots all features.
    #' @return A [ggplot2][ggplot2::ggplot] object
    plot_convergence = function(features = NULL) {
      require_package("ggplot2")

      if (is.null(self$convergence_history)) {
        cli::cli_abort("No convergence history available. Run $compute() first.")
      }

      # Create a copy to avoid modifying the original
      plot_data = copy(self$convergence_history)

      if (!is.null(features)) {
        plot_data = plot_data[feature %in% features]
      }

      p = ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = n_permutations, y = importance, fill = feature, color = feature)
      ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = importance - se, ymax = importance + se),
          alpha = 1 / 3
        ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(
          title = "SAGE Value Convergence",
          subtitle = if (self$converged) {
            sprintf(
              "Converged after %d permutations (saved %d)",
              self$n_permutations_used,
              self$n_permutations - self$n_permutations_used
            )
          } else {
            sprintf("Completed all %d permutations", self$n_permutations)
          },
          x = "Number of Permutations",
          y = "SAGE Value",
          color = "Feature",
          fill = "Feature"
        ) +
        ggplot2::theme_minimal(base_size = 14)

      if (self$converged) {
        p = p +
          ggplot2::geom_vline(
            xintercept = self$n_permutations_used,
            linetype = "dashed",
            color = "red",
            alpha = 0.5
          )
      }

      p
    }
  ),

  private = list(
    # This function computes the SAGE values for a single resampling iteration.
    # It iterates through permutations of features, evaluates coalitions, and calculates marginal contributions.
    .compute_sage_scores = function(
      learner,
      test_dt,
      n_permutations,
      batch_size = NULL,
      early_stopping = FALSE,
      se_threshold = 0.01,
      min_permutations = 10L,
      check_interval = 1L
    ) {
      # Initialize numeric vectors to store marginal contributions and their squares for variance calculation.
      # We track both sum and sum of squares to calculate running variance and standard errors.
      sage_values = numeric(length(self$features)) # Sum of marginal contributions
      sage_values_sq = numeric(length(self$features)) # Sum of squared marginal contributions
      names(sage_values) = self$features
      names(sage_values_sq) = self$features

      # Pre-generate `n_permutations` permutations upfront
      # Relevant for reproducibility, especially when using early stopping or parallel processing.
      # Example: if self$features = c("x1", "x2", "x3") and n_permutations = 2,
      # all_permutations might be list(c("x2", "x1", "x3"), c("x3", "x1", "x2"))
      all_permutations = replicate(n_permutations, sample(self$features), simplify = FALSE)

      # Initialize variables for iterative checkpoint-based computation.
      # This allows for early stopping based on convergence and provides progress updates.
      convergence_history = list() # Stores SAGE values at each checkpoint for convergence tracking
      n_completed = 0 # Number of permutations processed so far
      converged = FALSE # Flag to indicate if convergence has been detected
      baseline_loss = NULL # Loss of the empty coalition (model with no features / all features marginalized)

      # Calculate total checkpoints for progress tracking.
      # A checkpoint is a group of 'check_interval' permutations.
      total_checkpoints = ceiling(n_permutations / check_interval)
      current_checkpoint = 0

      # Start checkpoint-based progress bar if progress display is enabled.
      if (xplain_opt("progress")) {
        cli::cli_progress_bar(
          "Computing SAGE values",
          total = total_checkpoints
        )
      }

      # Main loop: Process permutations in checkpoints until all permutations are done or convergence is reached.
      while (n_completed < n_permutations && !converged) {
        # Determine the size of the current checkpoint.
        # This ensures that the last checkpoint processes only the remaining permutations.
        checkpoint_size = min(check_interval, n_permutations - n_completed)
        # Define the indices of permutations to be processed in this checkpoint.
        checkpoint_perms = (n_completed + 1):(n_completed + checkpoint_size)

        # Get the actual permutation sequences for this checkpoint from the pre-generated list.
        checkpoint_permutations = all_permutations[checkpoint_perms]

        # Build this checkpoint's growing-prefix coalitions. The
        # empty coalition is prepended only in the first checkpoint;
        # its loss is the baseline anchor for marginal contributions.
        # Same single-batch call/order as before, so the RNG-bearing
        # marginal sampling inside .evaluate_coalitions_batch is
        # byte-identical to the pre-refactor scheme.
        checkpoint_coalitions = sage_growing_coalitions(checkpoint_permutations)
        offset = 0L
        if (n_completed == 0) {
          checkpoint_coalitions = c(list(character(0)), checkpoint_coalitions)
          offset = 1L
        }

        # Progress: one tick per checkpoint (unchanged cadence).
        current_checkpoint = current_checkpoint + 1

        # Evaluate all coalitions collected in this checkpoint in a single batch.
        # This is a performance optimization to minimize prediction calls to the learner.
        checkpoint_losses = private$.evaluate_coalitions_batch(
          learner,
          test_dt,
          checkpoint_coalitions,
          batch_size
        )

        # Update progress bar.
        if (xplain_opt("progress")) {
          cli::cli_progress_update(inc = 1)
        }

        # Store the baseline loss (loss of the empty coalition) from the first checkpoint.
        # This is the model's performance when no features are available.
        if (n_completed == 0) {
          baseline_loss = checkpoint_losses[1] # The first element is always the empty coalition's loss
        }

        # Closed-form accumulation over the growing-prefix losses.
        # Every permutation is a full feature permutation, so each
        # coalition's loss index is computed directly; `offset` skips
        # the leading empty-coalition slot present in the first
        # checkpoint. Replaces the former O(n^2) which(sapply())
        # coalition-map lookup.
        acc = sage_marginal_contributions(
          checkpoint_permutations,
          checkpoint_losses,
          baseline_loss,
          self$features,
          offset = offset
        )
        # Name-aligned add (defensive: positional add is only valid if
        # orders match; reindex by name to be safe).
        sage_values = sage_values + acc$sv[names(sage_values)]
        sage_values_sq = sage_values_sq + acc$sv_sq[names(sage_values_sq)]

        # Update the count of completed permutations.
        n_completed = n_completed + checkpoint_size

        # Calculate the current average SAGE values and standard errors based on completed permutations.
        current_avg = sage_values / n_completed

        # Calculate running variance and standard errors for each feature
        # Variance = E[X^2] - E[X]^2, SE = sqrt(Var / n)
        current_variance = (sage_values_sq / n_completed) - (current_avg^2)
        # Ensure variance is non-negative (numerical precision issues)
        current_variance[current_variance < 0] = 0
        current_se = sqrt(current_variance / n_completed)

        if (xplain_opt("debug")) {
          cli::cli_alert_info("SAGE values after {.val {n_completed}} permutations")
          cli::cli_ol(c(
            "SAGE values: {.val {round(current_avg, 4)}}",
            "current SE: {.val {round(current_se, 3)}}",
            "Completed: {.val {n_completed}}"
          ))
        }

        # Store the current average SAGE values and standard errors in the convergence history.
        # Used for plotting, early stopping, and uncertainty quantification.
        checkpoint_history = data.table(
          n_permutations = n_completed,
          feature = names(current_avg),
          importance = as.numeric(current_avg),
          se = as.numeric(current_se)
        )
        convergence_history[[length(convergence_history) + 1]] = checkpoint_history

        # Check for convergence if early stopping is enabled and enough permutations have been processed.
        if (early_stopping && n_completed >= min_permutations && length(convergence_history) > 1) {
          # Get SAGE values from the current checkpoint.
          curr_checkpoint = convergence_history[[length(convergence_history)]]

          # Ensure features are in the same order for comparison.
          curr_checkpoint_ordered = copy(curr_checkpoint)[order(feature)]
          curr_importance_values = curr_checkpoint_ordered$importance
          curr_se_values = curr_checkpoint_ordered$se

          # Calculate range of importance values across all features (matching fippy)
          importance_range = max(curr_importance_values, na.rm = TRUE) -
            min(curr_importance_values, na.rm = TRUE)

          # Normalize SE by range to get relative SE (matching fippy's formula)
          # fippy: ratio = SE / range, convergence if max(ratio) < threshold
          # https://github.com/gcskoenig/fippy/blob/a7a37aa5511f7074ead3289c89b1ae80036982cb/src/fippy/explainers/utils.py#L40-L42
          if (importance_range > 0 && is.finite(importance_range)) {
            relative_se_values = curr_se_values / importance_range
            max_relative_se = max(relative_se_values, na.rm = TRUE)
          } else {
            # If range is 0 or invalid, use absolute SE (fallback)
            max_relative_se = max(curr_se_values, na.rm = TRUE)
          }

          # Check convergence: max relative SE below threshold
          converged = is.finite(max_relative_se) && max_relative_se < se_threshold
          convergence_msg = c(
            "v" = "SAGE converged after {.val {n_completed}} permutations",
            "i" = "Maximum relative SE: {.val {round(max_relative_se, 4)}} (threshold: {.val {se_threshold}})",
            "i" = "Saved {.val {n_permutations - n_completed}} permutations"
          )

          if (xplain_opt("verbose") && converged) {
            cli::cli_inform(convergence_msg)
          }
        }
      }

      # Close the progress bar.
      if (xplain_opt("progress")) {
        cli::cli_progress_done()
      }

      # Calculate the final average SAGE values based on all completed permutations.
      final_sage_values = sage_values / n_completed

      # Return the computed scores and convergence data.
      list(
        scores = data.table(
          feature = names(final_sage_values),
          importance = as.numeric(final_sage_values)
        ),
        convergence_data = list(
          convergence_history = if (length(convergence_history) > 0) {
            rbindlist(convergence_history)
          } else {
            NULL
          },
          converged = converged,
          n_permutations_used = n_completed
        )
      )
    },

    # Template method: Defines the complete prediction and aggregation pipeline
    # Subclasses only need to implement .expand_coalitions_data()
    .evaluate_coalitions_batch = function(learner, test_dt, all_coalitions, batch_size = NULL) {
      n_test = nrow(test_dt)

      if (xplain_opt("debug")) {
        cli::cli_inform("Evaluating {.val {length(all_coalitions)}} coalitions")
      }

      # STEP 1: Subclass-specific data expansion (abstract method)
      # combined data has rows `n_samples * nrow(test_dt) * length(all_coalitions)`
      # Full coalition -> return is just test_dt
      combined_data = private$.expand_coalitions_data(test_dt, all_coalitions)
      # STEPS 2-5: Shared processing pipeline using general utilities
      predictions = sage_batch_predict(
        learner,
        combined_data,
        self$task,
        batch_size,
        self$task$task_type
      )
      if (anyNA(predictions)) {
        cli::cli_warn("Encountered missing values in model prediction")
      }
      avg_preds = sage_aggregate_predictions(
        combined_data,
        predictions,
        self$task$task_type,
        self$task$class_names
      )

      # Private method (needs self$task and self$measure)
      coalition_losses = private$.calculate_coalition_losses(avg_preds, n_test, test_dt)

      # SAGE values are score reductions loss(empty) - loss(S). Negating the scores of a
      # measure that is maximized (e.g. classif.acc) keeps "positive = helps" for all measures.
      if (isFALSE(self$measure$minimize)) {
        coalition_losses = -coalition_losses
      }
      coalition_losses
    },

    # Abstract method - must be implemented by subclasses
    # Returns: data.table with all feature columns plus .coalition_id and .test_instance_id
    .expand_coalitions_data = function(test_dt, all_coalitions) {
      cli::cli_abort(c(
        "Abstract method not implemented",
        "i" = "Subclasses must implement {.fn .expand_coalitions_data}",
        "i" = "This method should return a data.table with:",
        "*" = "All feature columns (with marginalized features replaced/sampled)",
        "*" = "{.field .coalition_id}: integer identifying which coalition",
        "*" = "{.field .test_instance_id}: integer identifying original test instance"
      ))
    },

    # Private method - needs self$task and self$measure
    # Calculates losses from averaged predictions for each coalition.
    # Returns losses ordered by ascending `.coalition_id`, matching the
    # order coalitions were built (the closed-form accumulation indexes
    # this positionally).
    .calculate_coalition_losses = function(avg_preds, n_test, test_dt) {
      .coalition_id = .test_instance_id = NULL # data.table NSE NOTE tax

      coalition_ids = sort(unique(avg_preds$.coalition_id))
      truth = test_dt[[self$task$target_names]]

      # Key once and subset by key (binary search) per coalition,
      # rather than scanning `avg_preds[.coalition_id == i]` in a loop.
      # The key includes .test_instance_id so each block is
      # ordered to align with `truth` (test instance 1..n_test).
      setkey(avg_preds, .coalition_id, .test_instance_id)
      measure = self$measure
      is_classif = self$task$task_type == "classif"
      class_names = self$task$class_names

      # For pointwise (`obs_loss`) measures whose predict_type matches
      # the data we hold (regression response; classification prob), we
      # can call the measure's own mlr3measures function (`$fun`)
      # directly on (truth, response/prob), skipping the per-coalition
      # Prediction object. This reuses the mlr3 measure
      # implementation and does the
      # correct per-measure aggregation. Anything else (classification
      # response measures that need a prob->class step, measures needing
      # task/model context, non-decomposable measures like AUC) takes
      # the canonical Prediction$score() path, so correctness holds for
      # every measure.
      direct = !is.null(measure$fun) &&
        "obs_loss" %in% measure$properties &&
        ((!is_classif && measure$predict_type == "response") ||
          (is_classif && measure$predict_type == "prob"))

      score_block = if (direct && is_classif) {
        function(block) {
          measure$fun(truth = truth, prob = as.matrix(block[, .SD, .SDcols = class_names]))
        }
      } else if (direct) {
        function(block) measure$fun(truth = truth, response = block$avg_pred)
      } else if (is_classif) {
        function(block) {
          PredictionClassif$new(
            row_ids = seq_len(n_test),
            truth = truth,
            prob = as.matrix(block[, .SD, .SDcols = class_names])
          )$score(measure)
        }
      } else {
        function(block) {
          PredictionRegr$new(
            row_ids = seq_len(n_test),
            truth = truth,
            response = block$avg_pred
          )$score(measure)
        }
      }

      coalition_losses = numeric(length(coalition_ids))
      for (k in seq_along(coalition_ids)) {
        coalition_losses[k] = score_block(avg_preds[list(coalition_ids[k])])
      }

      coalition_losses
    }
  )
)
