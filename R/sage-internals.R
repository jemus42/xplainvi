#' Batch Predict for SAGE
#'
#' Performs batched prediction on combined data to manage memory usage.
#' Supports both classification (probability predictions) and regression.
#'
#' @param learner ([`Learner`][mlr3::Learner]) Trained mlr3 learner.
#' @param combined_data (`data.table`) Data with feature columns to predict on.
#' @param task ([`Task`][mlr3::Task]) mlr3 task object.
#' @param batch_size (`integer(1)` or `NULL`) Batch size for predictions. If `NULL` or if
#'   `total_rows <= batch_size`, processes all data at once.
#' @param task_type (`character(1)`) Task type, either `"classif"` or `"regr"`.
#'
#' @return For classification: `matrix` of class probabilities (n_rows x n_classes).
#'   For regression: `numeric` vector of predictions (length n_rows).
#'
#' @keywords internal
sage_batch_predict = function(learner, combined_data, task, batch_size, task_type) {
  total_rows = nrow(combined_data)

  if (!is.null(batch_size) && total_rows > batch_size) {
    # Batched prediction
    n_batches = ceiling(total_rows / batch_size)
    all_predictions = vector("list", n_batches)

    for (batch_idx in seq_len(n_batches)) {
      start_row = (batch_idx - 1) * batch_size + 1
      end_row = min(batch_idx * batch_size, total_rows)
      batch_data = combined_data[start_row:end_row]

      if (xplain_opt("debug")) {
        cli::cli_inform(
          "Predicting on {.val {nrow(batch_data)}} instances in batch {.val {batch_idx}/{n_batches}}"
        )
      }

      pred_result = if (is.function(learner$predict_newdata_fast)) {
        learner$predict_newdata_fast(newdata = batch_data, task = task)
      } else {
        learner$predict_newdata(newdata = batch_data, task = task)
      }

      all_predictions[[batch_idx]] = if (task_type == "classif") {
        pred_result$prob
      } else {
        pred_result$response
      }
    }

    # Combine predictions from all batches
    if (task_type == "classif") {
      do.call(rbind, all_predictions)
    } else {
      do.call(c, all_predictions)
    }
  } else {
    # Single prediction without batching
    if (xplain_opt("debug")) {
      cli::cli_inform("Predicting on {.val {nrow(combined_data)}} instances at once")
    }

    pred_result = if (is.function(learner$predict_newdata_fast)) {
      learner$predict_newdata_fast(newdata = combined_data, task = task)
    } else {
      learner$predict_newdata(newdata = combined_data, task = task)
    }

    if (task_type == "classif") {
      pred_result$prob
    } else {
      pred_result$response
    }
  }
}

#' Aggregate Predictions by Coalition and Test Instance
#'
#' Averages predictions across multiple samples (reference data or conditional samples)
#' for each unique combination of coalition and test instance.
#'
#' @param combined_data (`data.table`) Data with columns `.coalition_id`, `.test_instance_id`,
#'   and feature columns.
#' @param predictions (`matrix` or `numeric`) For classification: matrix of class probabilities.
#'   For regression: numeric vector of predictions.
#' @param task_type (`character(1)`) Task type, either `"classif"` or `"regr"`.
#' @param class_names (`character()` or `NULL`: `NULL`) Character vector of class names. Required
#'   for classification, ignored for regression.
#'
#' @return `data.table` with columns:
#'   - `.coalition_id`: Coalition identifier (integer)
#'   - `.test_instance_id`: Test instance identifier (integer)
#'   - For classification: One column per class with averaged probabilities (numeric)
#'   - For regression: `avg_pred` column with averaged predictions (numeric)
#'
#' @keywords internal
sage_aggregate_predictions = function(combined_data, predictions, task_type, class_names = NULL) {
  if (task_type == "classif") {
    # Add prediction columns to combined_data
    n_classes = ncol(predictions)
    for (j in seq_len(n_classes)) {
      combined_data[, paste0(".pred_class_", j) := predictions[, j]]
    }

    # Aggregate: calculate mean probability for each class, grouped by coalition and test instance
    agg_cols = paste0(".pred_class_", seq_len(n_classes))
    avg_preds = combined_data[,
      lapply(.SD, function(x) mean(x, na.rm = TRUE)),
      .SDcols = agg_cols,
      by = c(".coalition_id", ".test_instance_id")
    ]

    # Rename aggregated columns to original class names
    setnames(avg_preds, agg_cols, class_names)
    avg_preds
  } else if (task_type == "regr") {
    # Regression: add predictions and aggregate
    .prediction = NULL # the data.table NSE NOTE tax
    combined_data[, .prediction := predictions]

    combined_data[,
      list(avg_pred = mean(.prediction, na.rm = TRUE)),
      by = c(".coalition_id", ".test_instance_id")
    ]
  }
}

#' Build row-major growing-prefix coalitions for a permutation list
#'
#' For each permutation, emit its growing prefixes
#' (`perm[1]`, `perm[1:2]`, ..., `perm`). Row-major over
#' `(permutation, step)`. Pure; no evaluation, no RNG.
#'
#' @param perm_sublist (`list`) Feature-name permutations.
#' @return `list` of character vectors (coalitions).
#' @keywords internal
#' @noRd
sage_growing_coalitions = function(perm_sublist) {
  # Pre-allocate to the exact coalition count (sum of permutation lengths)
  coalitions = vector("list", sum(lengths(perm_sublist)))
  k = 1L
  for (i in seq_along(perm_sublist)) {
    perm = perm_sublist[[i]]
    for (j in seq_along(perm)) {
      coalitions[[k]] = perm[seq_len(j)]
      k = k + 1L
    }
  }
  coalitions
}

#' Accumulate SAGE marginal contributions from a loss vector
#'
#' Given growing-prefix losses laid out row-major over
#' `(permutation, step)` (optionally preceded by `offset` leading
#' slots, e.g. an empty-coalition entry), accumulate per-feature SAGE
#' value sums and squared sums. Every permutation is a full feature
#' permutation, so the loss index is closed-form: no search/map.
#' Pure; order-independent across permutations.
#'
#' @param perm_sublist (`list`) Feature-name permutations.
#' @param losses (`numeric`) Losses for `offset` leading slots then the
#'   row-major growing-prefix coalitions of `perm_sublist`.
#' @param baseline (`numeric(1)`) Empty-coalition loss anchor.
#' @param feature_names (`character`) Names for the output vectors.
#' @param offset (`integer(1)`: `0L`) Leading loss slots to skip.
#' @return `list(sv, sv_sq)` named numeric vectors over `feature_names`.
#' @keywords internal
#' @noRd
sage_marginal_contributions = function(perm_sublist, losses, baseline, feature_names, offset = 0L) {
  sv = numeric(length(feature_names))
  sv_sq = numeric(length(feature_names))
  names(sv) = feature_names
  names(sv_sq) = feature_names

  for (i in seq_along(perm_sublist)) {
    perm = perm_sublist[[i]]
    p = length(perm)
    prev_loss = baseline
    for (j in seq_len(p)) {
      feature = perm[j]
      current_loss = losses[offset + (i - 1L) * p + j]
      contribution = prev_loss - current_loss
      sv[feature] = sv[feature] + contribution
      sv_sq[feature] = sv_sq[feature] + contribution^2
      prev_loss = current_loss
    }
  }

  list(sv = sv, sv_sq = sv_sq)
}

# Number of evaluated coalitions: one empty-coalition baseline plus m growing prefixes per
# permutation, or all 2^m for exact enumeration. The currency in which estimators are
# comparable, unlike their own budget units.
sage_n_evals = function(estimator, m, budget) {
  if (identical(estimator, "exact")) 2^m else 1 + budget * m
}

sage_assert_exact_budget = function(m, max_features) {
  if (m > max_features) {
    cli::cli_abort(c(
      "The exact estimator would enumerate {.val {2^m}} coalitions of {m} features.",
      "i" = "This exceeds the {.arg max_features} cap ({max_features}); the cost grows as 2^n_features.",
      "i" = "Increase {.arg max_features} to override, or use {.code estimator = \"permutation\"} instead."
    ))
  }
  invisible(NULL)
}

# Point out a sampling budget that costs at least as much as enumeration. Skipped under
# early stopping (the budget is then an upper bound) and for very small or very large
# feature sets, where the comparison is moot.
sage_inform_budget_vs_exact = function(m, n_permutations, early_stopping = FALSE) {
  if (!xplain_opt("verbose") || early_stopping || m < 3L || m > 30L) {
    return(invisible(NULL))
  }
  n_exact = 2^m
  evals = sage_n_evals("permutation", m, n_permutations)
  if (evals >= n_exact) {
    cli::cli_inform(c(
      "i" = "The permutation estimator will evaluate {evals} coalitions,
             at least as many as enumerating all {n_exact} (2^{m}) coalitions.",
      "i" = "{.code estimator = \"exact\"} computes SAGE values without coalition-sampling
             error at the same or lower cost (for {.cls ConditionalSAGE}, oversampling can
             still be deliberate; see the {.arg estimator} docs)."
    ))
  }
  invisible(NULL)
}

# Convergence criterion of the reference Python `sage` package (`detect_convergence`):
# largest SE relative to the spread of the SAGE values. NA (never converged) if any SE is
# missing, e.g. before the second permutation.
sage_convergence_ratio = function(importance, se) {
  if (anyNA(importance) || anyNA(se)) {
    return(NA_real_)
  }
  spread = max(importance) - min(importance)
  # A degenerate spread (single feature, or all features equal) leaves nothing to normalize
  # by, so the absolute standard error is used instead.
  ratio = if (spread > 0 && is.finite(spread)) max(se) / spread else max(se)
  if (is.finite(ratio)) ratio else NA_real_
}
