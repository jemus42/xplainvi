# =============================================================================
# MarginalSAGE Tests
# =============================================================================

# -----------------------------------------------------------------------------
# Basic functionality
# -----------------------------------------------------------------------------

test_that("MarginalSAGE default behavior with minimal parameters", {
  # Use small params for test speed
  test_default_behavior(MarginalSAGE, task_type = "regr", n_permutations = 2L, n_samples = 20L)
})

test_that("MarginalSAGE works with classification tasks", {
  # Binary classification
  task_binary = tgen("2dnormals")$generate(n = 100)
  sage_binary = MarginalSAGE$new(
    task = task_binary,
    learner = lrn("classif.rpart", predict_type = "prob"),
    n_permutations = 2L
  )
  checkmate::expect_r6(sage_binary, c("FeatureImportanceMethod", "SAGE", "MarginalSAGE"))
  sage_binary$compute()
  expect_importance_dt(sage_binary$importance(), features = sage_binary$features)

  # Multiclass classification
  task_multi = tgen("cassini")$generate(n = 100)
  sage_multi = MarginalSAGE$new(
    task = task_multi,
    learner = lrn("classif.rpart", predict_type = "prob"),
    n_permutations = 2L
  )
  sage_multi$compute()
  expect_importance_dt(sage_multi$importance(), features = sage_multi$features)
  expect_length(task_multi$class_names, 3L)
})

test_that("MarginalSAGE featureless learner produces zero importance", {
  # Use small params for test speed
  test_featureless_zero_importance(
    MarginalSAGE,
    task_type = "regr",
    n_permutations = 2L,
    n_samples = 20L
  )
})

# -----------------------------------------------------------------------------
# Sensible results
# -----------------------------------------------------------------------------

test_that("MarginalSAGE friedman1 produces sensible ranking", {
  # Use small params for test speed
  test_friedman1_sensible_ranking(MarginalSAGE, n = 200L, n_permutations = 2L, n_samples = 20L)
})

# -----------------------------------------------------------------------------
# Resampling
# -----------------------------------------------------------------------------

test_that("MarginalSAGE with cross-validation resampling", {
  task = tgen("friedman1")$generate(n = 200)

  sage = MarginalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    resampling = rsmp("cv", folds = 3),
    n_permutations = 2L
  )
  sage$compute()

  expect_importance_dt(sage$importance(), features = sage$features)
  checkmate::expect_data_table(
    sage$scores(),
    types = c("integer", "character", "numeric"),
    nrows = sage$resampling$iters * length(sage$features),
    ncols = 3,
    any.missing = FALSE
  )
})

# -----------------------------------------------------------------------------
# Single feature
# -----------------------------------------------------------------------------

test_that("MarginalSAGE with single feature", {
  task = tgen("friedman1")$generate(n = 100)

  sage = MarginalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    features = "important4",
    n_permutations = 2L
  )
  sage$compute()

  expect_importance_dt(sage$importance(), features = "important4")
  expect_equal(nrow(sage$importance()), 1L)
})

# -----------------------------------------------------------------------------
# n_samples parameter
# -----------------------------------------------------------------------------

test_that("MarginalSAGE with custom n_samples", {
  task = tgen("friedman1")$generate(n = 200)

  sage = MarginalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    n_samples = 30L,
    n_permutations = 2L
  )
  sage$compute()

  expect_importance_dt(sage$importance(), features = sage$features)
})

# -----------------------------------------------------------------------------
# Reproducibility
# -----------------------------------------------------------------------------

test_that("MarginalSAGE reproducibility with same seed", {
  task = tgen("2dnormals")$generate(n = 100)
  learner = lrn("classif.rpart", predict_type = "prob")
  measure = msr("classif.ce")

  set.seed(42)
  sage1 = MarginalSAGE$new(
    task = task,
    learner = learner,
    measure = measure,
    n_permutations = 3L
  )
  sage1$compute()
  result1 = sage1$importance()

  set.seed(42)
  sage2 = MarginalSAGE$new(
    task = task,
    learner = learner,
    measure = measure,
    n_permutations = 3L
  )
  sage2$compute()
  result2 = sage2$importance()

  # Results should be identical with same seed
  expect_equal(result1$importance, result2$importance, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Parameter validation
# -----------------------------------------------------------------------------

test_that("MarginalSAGE parameter validation", {
  task = tgen("friedman1")$generate(n = 50)
  learner = lrn("regr.rpart")

  # n_permutations must be positive integer
  expect_error(MarginalSAGE$new(task = task, learner = learner, n_permutations = 0L))
  expect_error(MarginalSAGE$new(task = task, learner = learner, n_permutations = -1L))
})

test_that("MarginalSAGE requires predict_type='prob' for classification", {
  task = tgen("2dnormals")$generate(n = 50)

  # Should error for classification without predict_type = "prob"
  expect_error(
    MarginalSAGE$new(
      task = task,
      learner = lrn("classif.rpart", predict_type = "response")
    ),
    "Classification learners require probability predictions for SAGE."
  )
})

# -----------------------------------------------------------------------------
# Convergence tracking
# -----------------------------------------------------------------------------

test_that("MarginalSAGE SE tracking in convergence_history", {
  task = tgen("friedman1")$generate(n = 30)
  learner = lrn("regr.rpart")
  measure = msr("regr.mse")

  sage = MarginalSAGE$new(
    task = task,
    learner = learner,
    measure = measure,
    n_permutations = 6L,
    n_samples = 20L
  )

  # Compute with early stopping to get convergence history
  sage$compute(early_stopping = TRUE, se_threshold = 0.05, check_interval = 2L)

  # Check that convergence_history exists and has SE column
  expect_false(is.null(sage$convergence_history))
  expect_contains(colnames(sage$convergence_history), "se")

  # Check structure of convergence_history
  expected_cols = c("n_permutations", "feature", "importance", "se")
  expect_setequal(colnames(sage$convergence_history), expected_cols)

  # SE values should be non-negative and finite
  se_values = sage$convergence_history$se
  checkmate::expect_numeric(se_values, lower = 0, finite = TRUE)

  # For each feature, SE should be in a reasonable range
  for (feat in unique(sage$convergence_history$feature)) {
    feat_data = sage$convergence_history[feature == feat]
    feat_data = feat_data[order(n_permutations)]

    if (nrow(feat_data) > 1) {
      # Just check that SE values are in a reasonable range and not exploding
      expect_lt(max(feat_data$se), 10)
      expect_lt(max(abs(diff(feat_data$se))), 5)
    }
  }

  # All features should be represented in convergence history
  expect_setequal(
    unique(sage$convergence_history$feature),
    sage$features
  )
})

test_that("MarginalSAGE SE-based convergence detection", {
  skip_on_cran() # ~1s - tests early stopping feature, not core SAGE

  task = tgen("friedman1")$generate(n = 100)
  learner = lrn("regr.rpart")
  measure = msr("regr.mse")

  sage = MarginalSAGE$new(
    task = task,
    learner = learner,
    measure = measure,
    n_permutations = 10L,
    n_samples = 20L
  )

  # Test with very loose SE threshold (should trigger convergence easily)
  sage$compute(
    early_stopping = TRUE,
    se_threshold = 100.0,
    min_permutations = 5L,
    check_interval = 1L
  )

  # Should converge early because SE will be well below 100.0
  expect_true(sage$converged)
  expect_lte(sage$n_permutations_used, 10L)

  # Reset for next test
  sage$reset()

  # Test with very strict SE threshold (should not converge)
  sage$compute(
    early_stopping = TRUE,
    se_threshold = 0.001,
    min_permutations = 5L,
    check_interval = 1L
  )

  # With very strict SE threshold, should not converge early
  expect_false(sage$converged)

  # Test with moderate SE threshold
  sage$reset()

  sage$compute(
    early_stopping = TRUE,
    se_threshold = 0.1,
    min_permutations = 5L,
    check_interval = 1L
  )

  # Should have convergence history with SE tracking regardless of convergence
  expect_false(is.null(sage$convergence_history))
  expect_contains(colnames(sage$convergence_history), "se")
})

test_that("coalition-loss direct measure$fun path matches per-coalition $score() (regression)", {
  # Guards the direct measure$fun path in .calculate_coalition_losses
  # (used for obs_loss regression measures) against the canonical
  # per-coalition Prediction + $score() path.
  task = sim_dgp_independent(40)
  # learner is irrelevant here: the private method under test only uses
  # task + measure, so use the dependency-free featureless learner.
  learner = lrn("regr.featureless")

  n_test = 12L
  n_coal = 5L
  set.seed(606)
  truth = rnorm(n_test)
  test_dt = data.table(y = truth)
  avg_preds = CJ(.coalition_id = seq_len(n_coal), .test_instance_id = seq_len(n_test))
  avg_preds[, avg_pred := rnorm(.N)]

  for (meas_id in c("regr.mse", "regr.mae")) {
    sage = MarginalSAGE$new(task, learner, measure = msr(meas_id), n_permutations = 2L)
    priv = sage$.__enclos_env__$private

    # copy(): the method keys avg_preds by reference
    fast = priv$.calculate_coalition_losses(copy(avg_preds), n_test, test_dt)

    slow = vapply(
      seq_len(n_coal),
      function(i) {
        cd = avg_preds[.coalition_id == i][order(.test_instance_id)]
        PredictionRegr$new(
          row_ids = seq_len(n_test),
          truth = truth,
          response = cd$avg_pred
        )$score(msr(meas_id))
      },
      numeric(1)
    )

    expect_equal(fast, slow, info = meas_id)
  }
})

test_that("coalition-loss direct measure$fun path matches per-coalition $score() (classif prob)", {
  # classif.logloss is a prob obs_loss measure -> direct measure$fun path.
  task = tgen("2dnormals")$generate(40)
  learner = lrn("classif.featureless", predict_type = "prob")
  classes = task$class_names

  n_test = 12L
  n_coal = 4L
  set.seed(707)
  truth = factor(sample(classes, n_test, replace = TRUE), levels = classes)
  test_dt = data.table(class = truth)
  setnames(test_dt, "class", task$target_names)

  # synthetic averaged class probabilities per (coalition, test instance)
  avg_preds = CJ(.coalition_id = seq_len(n_coal), .test_instance_id = seq_len(n_test))
  pm = matrix(runif(nrow(avg_preds) * length(classes)), ncol = length(classes))
  pm = pm / rowSums(pm)
  for (j in seq_along(classes)) {
    avg_preds[, (classes[j]) := pm[, j]]
  }

  sage = MarginalSAGE$new(task, learner, measure = msr("classif.logloss"), n_permutations = 2L)
  priv = sage$.__enclos_env__$private
  fast = priv$.calculate_coalition_losses(copy(avg_preds), n_test, test_dt)

  slow = vapply(
    seq_len(n_coal),
    function(i) {
      cd = avg_preds[.coalition_id == i][order(.test_instance_id)]
      PredictionClassif$new(
        row_ids = seq_len(n_test),
        truth = truth,
        prob = as.matrix(cd[, .SD, .SDcols = classes])
      )$score(msr("classif.logloss"))
    },
    numeric(1)
  )

  expect_equal(fast, slow)
})

test_that("maximized measures yield the same SAGE values as their minimized counterpart", {
  # classif.acc = 1 - classif.ce per coalition, so the score reductions must coincide
  # once the maximized measure's scores are negated internally.
  task = tgen("2dnormals")$generate(n = 100)
  learner = lrn("classif.rpart", predict_type = "prob")
  resampling = rsmp("holdout")$instantiate(task)

  set.seed(2311)
  ce = MarginalSAGE$new(task, learner, msr("classif.ce"), resampling, n_permutations = 5L, n_samples = 20L)
  ce$compute()
  set.seed(2311)
  acc = MarginalSAGE$new(task, learner, msr("classif.acc"), resampling, n_permutations = 5L, n_samples = 20L)
  acc$compute()

  expect_equal(acc$importance()$importance, ce$importance()$importance, tolerance = 1e-10)
  expect_gt(sum(acc$importance()$importance), 0)
})

test_that("standardize = TRUE does not mutate stored scores", {
  set.seed(6389)
  task = sim_dgp_independent(n = 150)
  sage = MarginalSAGE$new(task, lrn("regr.rpart"), n_permutations = 3L, n_samples = 20L)
  sage$compute()
  before = data.table::copy(sage$scores())
  sage$importance(standardize = TRUE)
  sage$importance(standardize = TRUE)
  expect_equal(sage$scores(), before)
})
