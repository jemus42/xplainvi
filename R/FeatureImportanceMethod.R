#' Feature Importance Method Class
#'
#' @export
FeatureImportanceMethod = R6Class(
  "FeatureImportanceMethod",
  public = list(
    #' @field label (`character(1)`) Method label.
    label = NA_character_,
    #' @field task ([mlr3::Task])
    task = NULL,
    #' @field learner ([mlr3::Learner])
    learner = NULL,
    #' @field measure ([mlr3::Measure])
    measure = NULL,
    #' @field resampling ([mlr3::Resampling]), instantiated upon construction.
    resampling = NULL,
    #' @field resample_result ([mlr3::ResampleResult]) of the original `learner` and `task`, used for baseline scores.
    resample_result = NULL,
    #' @field features (`character`: `NULL`) Features of interest. By default, importances will be computed for each feature
    #'   in `task`, but optionally this can be restricted to at least one feature. Ignored if `groups` is specified.
    features = NULL,
    #' @field groups (`list`: `NULL`) A (named) list of features (names or indices as in `task`).
    #'   If `groups` is specified, `features` is ignored.
    #'   Importances will be calculated for group of features at a time, e.g., in [PFI] not one but the group of features will be permuted at each step.
    #'   Analogously in [WVIM], each group of features will be left out (or in) for each model refit.
    #'   Not all methods support groups (e.g., [SAGE]).
    groups = NULL,
    #' @field param_set ([paradox::ps()])
    param_set = ps(),
    #' @field predictions ([data.table][data.table::data.table]) Feature-specific prediction objects provided for some methods ([PFI], [WVIM]). Contains columns for feature of interest, resampling iteration, refit or perturbation iteration, and [mlr3::Prediction] objects.
    predictions = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #' This is typically intended for use by derived classes.
    #' @param task,learner,measure,resampling,features,groups,param_set,label Used to set fields
    initialize = function(
      task,
      learner,
      measure = NULL,
      resampling = NULL,
      features = NULL,
      groups = NULL,
      param_set = paradox::ps(),
      label
    ) {
      self$task = mlr3::assert_task(task)
      self$learner = mlr3::assert_learner(learner, task = task, task_type = task$task_type)

      if (is.null(measure)) {
        self$measure = switch(
          task$task_type,
          "classif" = mlr3::msr("classif.ce"),
          "regr" = mlr3::msr("regr.mse")
        )
        if (xplain_opt("verbose")) {
          cli::cli_alert_info(
            "No {.cls Measure} provided, using {.code measure = msr(\"{self$measure$id}\")}"
          )
        }
      } else {
        self$measure = mlr3::assert_measure(measure, task = task, learner = learner)
      }
      self$param_set = paradox::assert_param_set(param_set)
      self$label = checkmate::assert_string(label, min.chars = 1)

      # Check features / groups
      # Default to using features, unless groups is specified
      if (is.null(groups)) {
        checkmate::assert_subset(features, self$task$feature_names, empty.ok = TRUE)
        self$features = features %||% self$task$feature_names
      } else {
        self$groups = check_groups(groups, all_features = self$task$feature_names)
        # check_groups ensures this produces a unique character vector
        self$features = unlist(groups, use.names = FALSE)
      }

      # resampling: default to holdout with default ratio if NULL
      if (is.null(resampling)) {
        resampling = mlr3::rsmp("holdout", ratio = 2 / 3)$instantiate(task)
        if (xplain_opt("verbose")) {
          cli::cli_inform(c(
            i = "No {.cls Resampling} provided, using {.code resampling = rsmp(\"holdout\", ratio = 2/3)} (test set size: {.val {length(resampling$test_set(1))}})"
          ))
        }
      } else {
        # Clone the resampling to avoid instantiating the resampling in the user's workspace
        resampling = mlr3::assert_resampling(resampling)$clone()

        # A pretrained learner requires a user-provided instantiated resampling
        # to define the test set explicitly. Auto-instantiation would pick an
        # arbitrary split unrelated to how the learner was trained.
        if (!is.null(learner$model) && !resampling$is_instantiated) {
          cli::cli_abort(c(
            "A pre-trained {.cls Learner} requires an instantiated {.cls Resampling}",
            i = "Instantiate the {.cls Resampling} before passing it, e.g. {.code rsmp(\"holdout\")$instantiate(task)}"
          ))
        }
      }
      if (!resampling$is_instantiated) {
        resampling$instantiate(task)
      }
      self$resampling = resampling

      # Check pretrained learner compatibility (multi-fold with pretrained learner)
      assert_pretrained(self$learner, self$task, self$resampling)
    },

    #' @description
    #' Compute feature importance scores
    #' @param store_backends (`logical(1): TRUE`) Whether to store backends.
    compute = function(store_backends = TRUE) {
      stop("Abstract method. Use a concrete implementation.")
    },

    #' @description
    #' Get aggregated importance scores.
    #' The stored [`measure`][mlr3::Measure] object's `aggregator` (default: `mean`) will be used to aggregated importance scores
    #' across resampling iterations and, depending on the method use, permutations ([PerturbationImportance] or refits [LOCO]).
    #' @param relation (character(1)) How to relate perturbed scores to originals ("difference" or "ratio").
    #'   If `NULL`, uses stored parameter value. This is only applicable for methods where importance is based on some
    #'   relation between baseline and post-modification loss, i.e. [PerturbationImportance] methods such as [PFI] or [WVIM] / [LOCO].
    #'   Not available for [SAGE] methods.
    #'   For `"ratio"`, features whose baseline score is `0` yield `NA` (with a warning) and therefore aggregate to `NA`.
    #' @param standardize (`logical(1)`: `FALSE`) If `TRUE`, importances are standardized by the highest score so all scores fall in `[-1, 1]`.
    #' @param ci_method (`character(1)`: `"none"`) Which confidence interval estimation method to use, defaulting to omitting
    #'   variance estimation (`"none"`).
    #'   If `"raw"`, uncorrected (too narrow) CIs are provided purely for informative purposes.
    #'   If `"nadeau_bengio"`, variance correction is performed according to Nadeau & Bengio (2003) as suggested by Molnar et al. (2023).
    #'   If `"quantile"`, empirical quantiles are used to construct confidence-like intervals.
    #'   These methods are model-agnostic and rely on suitable `resampling`s, e.g. subsampling with 15 repeats for `"nadeau_bengio"`.
    #'   See details.
    #' @param conf_level (`numeric(1)`: `0.95`) Confidence level to use for confidence interval construction when `ci_method != "none"`.
    #' @param alternative (`character(1)`: `"two.sided"`) Type of alternative hypothesis for statistical tests.
    #'   `"greater"` tests H0: importance <= 0 vs H1: importance > 0 (one-sided).
    #'   `"two.sided"` tests H0: importance = 0 vs H1: importance != 0.
    #'   Only used when `ci_method != "none"`.
    #' @param ... Additional arguments passed to specialized methods, if any.
    #' @return ([data.table][data.table::data.table]) Aggregated importance scores with columns `"feature"`, `"importance"`,
    #' and depending on `ci_method` also `"se"`, `"statistic"`, `"p.value"`, `"conf_lower"`, `"conf_upper"`.
    #'
    #' @details
    #'
    #' ## Confidence Interval Methods
    #'
    #' The parametric methods (`"raw"`, `"nadeau_bengio"`) return standard error (`se`),
    #' test statistic (`statistic`), p-value (`p.value`), and confidence bounds
    #' (`conf_lower`, `conf_upper`). The `"quantile"` method returns only lower and upper bounds.
    #'
    #' **`"raw"`: Uncorrected (!) t-test**
    #' Uses a standard t-test assuming independence of resampling iterations.
    #' - SE = sd(resampling scores) / sqrt(n_iters)
    #' - Test statistic: t = importance / SE with df = n_iters - 1
    #' - P-value: From t-distribution (one-sided or two-sided depending on `alternative`)
    #' - CIs: importance +/- qt(1 - alpha, df) * SE
    #'
    #' **Warning**: These CIs are too narrow because resampling iterations share
    #' training data and are not independent.
    #' This method is included only for demonstration purposes.
    #'
    #' **`"nadeau_bengio"`: Corrected t-test**
    #' Applies the Nadeau & Bengio (2003) correction to account for correlation between
    #' resampling iterations due to overlapping training sets.
    #' - Correction factor: (1/n_iters + n_test/n_train)
    #' - SE = sqrt(correction_factor * var(resampling scores))
    #' - Test statistic and p-value: As in `"raw"`, but with corrected SE
    #'
    #' Recommended with bootstrap or subsampling (>= 10 iterations).
    #'

    #' **`"quantile"`: Non-parametric empirical method**
    #' Uses the resampling distribution directly without parametric assumptions.
    #' - CIs: Empirical quantiles of the resampling distribution
    #'
    #' This method does not provide `se`, `statistic`, or `p.value`.
    #'
    #' ## Method-Specific CI Methods
    #'
    #' Some importance methods provide additional CI methods tailored to their approach:
    #'
    #' - **[CFI]**: Adds `"cpi"` (Conditional Predictive Impact), which uses observation-wise
    #'   loss differences with holdout resampling. Supports t-test, Wilcoxon, Fisher permutation,
    #'   and binomial tests. See Watson & Wright (2021).
    #'
    #' ## Practical Recommendations
    #'
    #' Variance estimates for importance scores are biased due to the resampling procedure.
    #' Molnar et al. (2023) suggest using the Nadeau & Bengio correction with approximately
    #' 15 iterations of subsampling.
    #'
    #' Bootstrapping can cause information leakage with learners that bootstrap internally
    #' (e.g., Random Forests), as observations may appear in both train and test sets.
    #' Prefer subsampling in such cases:
    #'
    #' ```r
    #' PFI$new(
    #'   task = sim_dgp_interactions(n = 1000),
    #'   learner = lrn("regr.ranger", num.trees = 100),
    #'   measure = msr("regr.mse"),
    #'   resampling = rsmp("subsampling", repeats = 15),
    #'   n_repeats = 20
    #' )
    #' ```
    #'
    #' The `"nadeau_bengio"` correction was validated for PFI; its use with other methods
    #' like LOCO or SAGE is experimental.
    #'
    #' @param p_adjust (`character(1)`: `"none"`) Method for p-value adjustment for multiple comparisons.
    #'   Accepts any method supported by [stats::p.adjust.methods], e.g. `"holm"`, `"bonferroni"`, `"BH"`, `"none"`.
    #'   Applied to p-values from `"raw"` and `"nadeau_bengio"` methods.
    #'   When `"bonferroni"`, confidence intervals are also adjusted (alpha/k).
    #'   For other correction methods (e.g. `"holm"`, `"BH"`), only p-values are adjusted;
    #'   confidence intervals remain at the nominal `conf_level` because these sequential/adaptive
    #'   procedures do not have a clean per-comparison alpha for CI construction.
    #'
    #' @references
    #' `r print_bib("nadeau_2003")`
    #' `r print_bib("molnar_2023")`
    #'
    importance = function(
      relation = NULL,
      standardize = FALSE,
      ci_method = c("none", "raw", "nadeau_bengio", "quantile"),
      conf_level = 0.95,
      alternative = c("two.sided", "greater"),
      p_adjust = "none",
      ...
    ) {
      if (is.null(private$.scores)) {
        cli::cli_inform(c(
          x = "No importances computed yet!"
        ))
        return(invisible(NULL))
      }

      # Catch unknown arguments that were not consumed by subclass methods
      if (...length() > 0) {
        dots = list(...)
        cli::cli_abort(c(
          "Unknown argument{?s}: {.arg {names(dots)}}.",
          i = "These arguments are not used by {.fun $importance} with {.code ci_method = \"{ci_method}\"}."
        ))
      }

      # Validate ci_method
      if (length(ci_method) > 1) {
        ci_method = ci_method[1]
      }
      checkmate::assert_choice(ci_method, choices = private$.ci_methods)
      checkmate::assert_number(conf_level, lower = 0, upper = 1)
      checkmate::assert_choice(p_adjust, choices = stats::p.adjust.methods)
      alternative = match.arg(alternative)

      # Get aggregator and scores
      aggregator = self$measure$aggregator %||% mean
      scores = self$scores(relation = relation)

      # Standardize first so variance calculations use standardized values
      if (standardize) {
        scores[, importance := importance / max(abs(importance), na.rm = TRUE)]
      }

      # Dispatch to appropriate aggregation function
      agg_importance = switch(
        ci_method,
        none = importance_none(scores, aggregator, conf_level),
        raw = importance_raw(
          scores,
          aggregator,
          conf_level,
          alternative,
          self$resample_result$iters,
          p_adjust = p_adjust
        ),
        nadeau_bengio = importance_nadeau_bengio(
          scores,
          aggregator,
          conf_level,
          alternative,
          self$resampling,
          self$resample_result$iters,
          p_adjust = p_adjust
        ),
        quantile = importance_quantile(scores, aggregator, conf_level, alternative),
        cli::cli_abort(c(
          "Variance method {.val {ci_method}} not found.",
          i = "Available methods: {.val {private$.ci_methods}}"
        ))
      )

      setkeyv(agg_importance, "feature")
      agg_importance[]
    },

    #' @description
    #' Calculate observation-wise importance scores.
    #'
    #' Requires that `$compute()` was run and that `measure` is decomposable and
    #' has an observation-wise loss (`Measure$obs_loss()`) associated with it.
    #' This is not the case for measure like `classif.auc`, which is not decomposable.
    #'
    #' @param relation (character(1)) How to relate perturbed scores to originals ("difference" or "ratio"). If `NULL`, uses stored parameter value. This is only applicable for methods where importance is based on some
    #' relation between baseline and post-modification loss, i.e. [PerturbationImportance] methods such as [PFI] or [WVIM] / [LOCO]. Not available for [SAGE] methods.
    #'   For `"ratio"`, features whose baseline score is `0` yield `NA` (with a warning) and therefore aggregate to `NA`.
    #'
    #' @return ([data.table][data.table::data.table]) Observation-wise losses and importance scores with columns
    #'   `"feature"`, `"iter_rsmp"`, `"iter_repeat"` (if applicable), `"row_ids"`, `"loss_baseline"`, `"loss_post"`, and `"obs_importance"`.
    obs_loss = function(relation = NULL) {
      if (!has_obs_loss(self$measure)) {
        cli::cli_warn(c(
          x = "{.cls Measure} {.val {self$measure$id}} does not have an observation-wise loss:",
          i = "Is it decomposable?"
        ))
        return(invisible(NULL))
      }
      if (is.null(private$.obs_losses)) {
        cli::cli_warn(c(
          x = "No observation-wise losses stored!",
          i = "Did you run {.fun $compute}?",
          i = "Not all methods support observation-wise losses"
        ))
        return(invisible(NULL))
      }

      relation = resolve_param(relation, self$param_set$values$relation, "difference")

      # Prepare baseline losses
      obs_loss_baseline = self$resample_result$obs_loss(measures = self$measure)
      # obs_loss_baseline[, let(truth = NULL, response = NULL)]
      setnames(
        obs_loss_baseline,
        old = c("iteration", self$measure$id),
        new = c("iter_rsmp", "loss_baseline")
      )

      obs_loss_combined = obs_loss_baseline[
        private$.obs_losses,
        on = .(iter_rsmp, row_ids),
        allow.cartesian = TRUE
      ]

      obs_loss_combined[,
        obs_importance := private$.compute_score(
          loss_baseline,
          loss_post,
          relation = relation
        )
      ]

      # Select / reorder column names, some may be
      # specific to methods and may not be present
      names_to_keep = c(
        "feature",
        "iter_rsmp",
        "iter_repeat",
        "row_ids",
        "loss_baseline",
        "loss_post",
        "obs_importance"
      )
      names_to_keep = intersect(names_to_keep, colnames(obs_loss_combined))

      obs_loss_combined[, .SD, .SDcols = names_to_keep][]
    },

    #' @description
    #' Resets all stored fields populated by `$compute`: `$resample_result`, `$scores`, `$obs_losses`, and `$predictions`.
    reset = function() {
      self$resample_result = NULL
      private$.scores = NULL
      private$.obs_losses = NULL
      self$predictions = NULL
      # SAGE-specific fields (only reset if they exist)
      if ("n_permutations_used" %in% names(self)) {
        self$n_permutations_used = NULL
      }
    },

    #' @description
    #' Print importance scores
    #'
    #' @param ... Passed to `print()`
    print = function(...) {
      cli::cli_h2(self$label)
      cli::cli_ul()
      cli::cli_li("Learner: {.val {self$learner$id}}")
      cli::cli_li("Task: {.val {self$task$id}}")
      if (is.null(self$groups)) {
        cli::cli_li(
          "{.emph {length(self$features)}} feature{?s} of interest: {.val {self$features}}"
        )
      } else {
        cli::cli_li("{.emph {length(self$groups)}} feature group{?s} of interest:")
        ol = cli::cli_ol()
        for (i in seq_along(groups)) {
          cli::cli_li("{.strong {names(groups)[i]}}: {.val {groups[[i]]}}")
        }
        cli::cli_end(ol)
      }
      cli::cli_li("Resampling: {.val {self$resampling$id}} ({.val {self$resampling$iters}} iters)")

      cli::cli_li("Parameters:")

      pv = self$param_set$values
      pidx = seq_along(pv)
      sapply(pidx, \(i) {
        cli::cli_ul("{.code {names(pv)[i]}}: {.val {pv[i]}}")
      })
      cli::cli_end()
      self$importance()
    },

    #' @description
    #' Calculate importance scores for each resampling iteration and sub-iterations
    #' (`iter_rsmp` in [PFI] for example).
    #'
    #' Iteration-wise importance are computed on the fly depending on the chosen relation
    #' (`difference` or `ratio`) to avoid re-computation if only a different relation is needed.
    #'
    #' @param relation (character(1)) How to relate perturbed scores to originals ("difference" or "ratio"). If `NULL`, uses stored parameter value. This is only applicable for methods where importance is based on some
    #' relation between baseline and post-modification loss, i.e. [PerturbationImportance] methods such as [PFI] or [WVIM] / [LOCO]. Not available for [SAGE] methods.
    #'   For `"ratio"`, features whose baseline score is `0` yield `NA` (with a warning) and therefore aggregate to `NA`.
    #'
    #' @return ([data.table][data.table::data.table]) Iteration-wise importance scores with columns for
    #'   `"feature"`, iteration indices, baseline and post-modification scores, and `"importance"`.
    scores = function(relation = NULL) {
      if (is.null(private$.scores)) {
        cli::cli_warn(c(
          x = "No importances computed yet!",
          i = "Did you run {.fun $compute}?"
        ))
        return(invisible(NULL))
      }
      if ("importance" %in% colnames(private$.scores)) {
        # If there is already an importance variable in the stored scores like in SAGE,
        # we can't calculate pre/post scores like in PFI, LOCO etc,
        # individual "scores" would have different meaning there.
        # Return a copy so callers (e.g. standardization in $importance()) cannot
        # mutate the stored scores by reference.
        return(data.table::copy(private$.scores))
      }

      relation = resolve_param(relation, self$param_set$values$relation, "difference")

      scores = data.table::copy(private$.scores)[,
        importance := private$.compute_score(
          score_baseline,
          score_post,
          relation = relation
        )
      ]

      setnames(
        scores,
        old = c("score_baseline", "score_post"),
        new = c(paste0(self$measure$id, c("_baseline", "_post")))
      )

      scores[]
    }
  ),
  private = list(
    # Registry of available variance methods
    .ci_methods = c("none", "raw", "nadeau_bengio", "quantile"),

    # Take the raw predictions as returned by $predict_newdata_fast and convert to Prediction object fitting the task type to simplify type-specific handling
    # @param raw_prediction `list` with elements `reponse` (vector) or `prob` (matrix) depending on task type.
    # @param test_row_ids `integer()` test set row ids, important to ensure predictions can be matched with original observations / baseline predictions
    # .construct_pred = function(raw_prediction, test_row_ids) {
    # 	truth = self$task$truth(rows = test_row_ids)

    # 	switch(
    # 		self$task$task_type,
    # 		classif = PredictionClassif$new(
    # 			row_ids = test_row_ids,
    # 			truth = truth,
    # 			response = raw_prediction$response, # vector of class names or NULL
    # 			prob = raw_prediction$prob # matrix for predict_type prob or NULL
    # 		),
    # 		regr = PredictionRegr$new(
    # 			row_ids = test_row_ids,
    # 			truth = truth,
    # 			response = raw_prediction$response # numeric
    # 		)
    # 	)
    # },

    # Utility to convert named list of groups of features into data.table to
    # make it a little easier to match group names and features in list columns etc
    # Used in WVIM where mlr3fselect stores "left in" features as list columns
    .groups_tbl = function() {
      group_tbl = data.table::data.table(
        group = names(self$groups),
        features_lst = unname(self$groups)
      )
      group_tbl[, features := vapply(features_lst, \(x) paste0(x, collapse = ";"), character(1))]
      group_tbl
    },

    # Scoring utility for computing importances
    #
    # Computes the `relation` of score before a change (e.g. PFI, LOCO, ...) and after.
    # If `minimize == TRUE`, then `scores_post - scores_pre` is computed for
    # `relation == "difference"`, otherwise `scores_pre - scores_post` is given.
    # If `minimize == FALSE`, then the order is flipped, ensuring that "higher value" means "more important".
    # @param scores_pre,scores_post (`numeric()`) Vector of scores or loss values at baseline / before (`_pre`) a modification, and after (`_post`) a modification (e.g., permutation or refit).
    # @param relation (`character(1)`: `"difference"`) Calculate the difference or `"ratio"` between pre and post modification value.
    .compute_score = function(
      scores_pre,
      scores_post,
      relation = c("difference", "ratio")
    ) {
      checkmate::assert_numeric(scores_pre, any.missing = FALSE)
      checkmate::assert_numeric(scores_post, any.missing = FALSE)
      checkmate::assert_true(length(scores_pre) == length(scores_post))
      relation = match.arg(relation)
      minimize = self$measure$minimize

      # General idea assuming a important feature:
      # For PFI and MSE (minimize == TRUE): post - baseline gives large value -> high importance
      # =="== and classif.acc (minimize == FALSE) -> post - baseline  = negative, so we flip
      # For WVIM when we "leave-in", the baseline scores are the empty model,
      # so we need to flip directions of the comparison as well to ensure "higher importance value" -> "more important"
      if (identical(self$direction, "leave-in")) {
        minimize = !minimize
      }

      # I know this could be more concise but for the time I prefer it to be very obvious in what happens when
      # General expectation -> higher score => more important
      if (minimize) {
        # Lower is better, e.g. ce, where scores_pre is expected to be smaller and scores_post larger
        switch(
          relation,
          difference = scores_post - scores_pre,
          ratio = safe_ratio(scores_post, scores_pre)
        )
      } else {
        # Higher is better, e.g. accuracy, where scores_pre is expected to be larger and scores_post smaller
        switch(
          relation,
          difference = scores_pre - scores_post,
          ratio = safe_ratio(scores_pre, scores_post)
        )
      }
    },

    # @field .scores ([data.table][data.table::data.table]) Iteration-wise importances scores. Essentially an aggregated form of .obs_losses (which may not be available), used as basis for the calculation in `$importance() `and `$scores()`.
    .scores = NULL,

    # @field .obs_losses ([data.table][data.table::data.table]) Observation-wise losses when available. Contains columns for row_ids, feature, iteration indices, individual loss values, and observation-wise losses for baseline and modified case.
    .obs_losses = NULL
  )
)
