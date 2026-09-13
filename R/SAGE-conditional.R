#' @title Conditional SAGE
#'
#' @description [SAGE] with conditional sampling (features are "marginalized" conditionally).
#' Uses [ConditionalARFSampler] as default [ConditionalSampler].
#'
#' @seealso [MarginalSAGE]
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE) && requireNamespace("mlr3learners", quietly = TRUE) && requireNamespace("arf", quietly = TRUE)
#' library(mlr3)
#' task = tgen("friedman1")$generate(n = 100)
#'
#' \donttest{
#' # Using default ConditionalARFSampler (also handles all mixed data)
#' sage = ConditionalSAGE$new(
#'   task = task,
#'   learner = lrn("regr.ranger", num.trees = 50),
#'   measure = msr("regr.mse"),
#'   n_permutations = 3L,
#'   n_samples = 20
#' )
#' sage$compute()
#' }
#' \donttest{
#' # For alternative conditional samplers:
#' custom_sampler = ConditionalGaussianSampler$new(
#'   task = task
#' )
#' sage_custom = ConditionalSAGE$new(
#'   task = task,
#'   learner = lrn("regr.ranger", num.trees = 50),
#'   measure = msr("regr.mse"),
#'   n_permutations = 5L,
#'   n_samples = 20,
#'   sampler = custom_sampler
#' )
#' sage_custom$compute()
#' }
#' @export
ConditionalSAGE = R6Class(
  "ConditionalSAGE",
  inherit = SAGE,
  public = list(
    #' @field sampler ([ConditionalSampler]) Sampler for conditional marginalization.
    sampler = NULL,

    #' @description
    #' Creates a new instance of the ConditionalSAGE class.
    #' @param task,learner,measure,resampling,features,n_permutations,batch_size,n_samples,early_stopping,se_threshold,min_permutations,check_interval Passed to [SAGE].
    #' @param sampler ([ConditionalSampler]) Optional custom sampler. Defaults to [ConditionalARFSampler].
    initialize = function(
      task,
      learner,
      measure = NULL,
      resampling = NULL,
      features = NULL,
      n_permutations = 10L,
      sampler = NULL,
      batch_size = 5000L,
      n_samples = 100L,
      early_stopping = FALSE,
      se_threshold = 0.025,
      min_permutations = 10L,
      check_interval = 1L
    ) {
      # Use ConditionalARFSampler by default
      if (is.null(sampler)) {
        sampler = ConditionalARFSampler$new(task)
        if (xplain_opt("verbose")) {
          cli::cli_alert_info(
            "No {.cls ConditionalSampler} provided, using {.cls ConditionalARFSampler} with default settings."
          )
        }
      } else {
        checkmate::assert_class(sampler, "ConditionalSampler")
      }

      super$initialize(
        task = task,
        learner = learner,
        measure = measure,
        resampling = resampling,
        features = features,
        n_permutations = n_permutations,
        batch_size = batch_size,
        n_samples = n_samples,
        early_stopping = early_stopping,
        se_threshold = se_threshold,
        min_permutations = min_permutations,
        check_interval = check_interval
      )

      # Store sampler (specific to ConditionalSAGE)
      self$sampler = sampler

      self$label = "Conditional SAGE"
    }
  ),

  private = list(
    # Conditional-specific data expansion: Sample from conditional distribution
    # For each coalition, sample marginalized features conditioned on coalition features
    # Uses n_samples to create multiple samples per test instance for averaging
    .expand_coalitions_data = function(test_dt, all_coalitions) {
      n_test = nrow(test_dt)
      all_coalition_data = vector("list", length(all_coalitions))

      # For each coalition, sample from conditional distribution P(X_{-S} | X_S)
      for (i in seq_along(all_coalitions)) {
        coalition = all_coalitions[[i]]
        marginalize_features = setdiff(self$features, coalition)

        if (length(marginalize_features) > 0) {
          n_samples = self$param_set$values$n_samples

          # Sample conditionally using a sampler. Pass `samples_per_row = n_samples`
          # on the unique test rows so the sampler can use its native efficient
          # multi-draw path (e.g. ConditionalARFSampler invokes
          # `arf::forge(n_synth = n_samples)` on n_test evidence rows rather than
          # n_samples * n_test replicated rows).
          # Returns test_dt with "marginalized" features replaced by conditional samples
          # sampler: P(marginalize_features | coalition_features)
          marginalized_test = self$sampler$sample_newdata(
            feature = marginalize_features,
            newdata = test_dt,
            conditioning_set = coalition,
            samples_per_row = n_samples
          )

          # Tag rows with their originating test instance. Sampler output is
          # draw-major: rows 1..n_test = draw 1 across all test rows in order,
          # rows n_test+1..2*n_test = draw 2, etc. Hence `times`, not `each`.
          # Must be done AFTER sampling since sampler does not preserve extra columns.
          marginalized_test[,
            .test_instance_id := rep(seq_len(n_test), times = n_samples)
          ]
        } else {
          # If marginalize_features is empty then we evaluate the coalition of all features
          # in the task, so there's no room for any sampling and we just take all data
          # no need to duplicate+average anything since prediction is deterministic (usually (I hope (right?)))
          marginalized_test = copy(test_dt)
          marginalized_test[, .test_instance_id := seq_len(n_test)]
        }

        # Add coalition ID for tracking which coalition this data belongs to
        marginalized_test[, .coalition_id := i]

        all_coalition_data[[i]] = marginalized_test
      }

      # Combine all coalitions into single data.table for batch prediction
      rbindlist(all_coalition_data)
    }
  )
)
