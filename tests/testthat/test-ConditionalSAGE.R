# =============================================================================
# ConditionalSAGE Tests
# =============================================================================
#
# Strategy for test speed:
# - Use ConditionalGaussianSampler (fast) for most tests
# - Use sim_dgp_* functions (4-5 features) instead of friedman1 (10 features)
# - Only use ARF sampler in tests specifically testing ARF behavior
# - ARF-specific tests get skip_on_cran()

# -----------------------------------------------------------------------------
# Basic functionality
# -----------------------------------------------------------------------------

test_that("ConditionalSAGE basic workflow with Gaussian sampler", {
  # Use small task (4 features) and fast Gaussian sampler
  task = sim_dgp_correlated(n = 100)
  sampler = ConditionalGaussianSampler$new(task)

  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    sampler = sampler,
    n_permutations = 2L,
    n_samples = 20L
  )

  checkmate::expect_r6(sage, c("FeatureImportanceMethod", "SAGE", "ConditionalSAGE"))
  sage$compute()
  expect_importance_dt(sage$importance(), features = sage$features)
})

test_that("ConditionalSAGE works with classification tasks", {
  # Binary classification - 2dnormals has 2 features
  task_binary = tgen("2dnormals")$generate(n = 50)
  sampler = ConditionalGaussianSampler$new(task_binary)

  sage_binary = ConditionalSAGE$new(
    task = task_binary,
    learner = lrn("classif.rpart", predict_type = "prob"),
    sampler = sampler,
    n_permutations = 2L,
    n_samples = 20L
  )
  checkmate::expect_r6(sage_binary, c("FeatureImportanceMethod", "SAGE", "ConditionalSAGE"))
  sage_binary$compute()
  expect_importance_dt(sage_binary$importance(), features = sage_binary$features)
})

test_that("ConditionalSAGE multiclass classification", {
  skip_on_cran() # multiclass with 3 features is slower
  task_multi = tgen("cassini")$generate(n = 50)
  sampler = ConditionalGaussianSampler$new(task_multi)

  sage_multi = ConditionalSAGE$new(
    task = task_multi,
    learner = lrn("classif.rpart", predict_type = "prob"),
    sampler = sampler,
    n_permutations = 2L,
    n_samples = 20L
  )
  sage_multi$compute()
  expect_importance_dt(sage_multi$importance(), features = sage_multi$features)
  expect_length(task_multi$class_names, 3L)
})

test_that("ConditionalSAGE featureless learner produces zero importance", {
  # Use small task (4 features) and fast Gaussian sampler
  task = sim_dgp_correlated(n = 100)
  sampler = ConditionalGaussianSampler$new(task)

  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.featureless"),
    measure = msr("regr.mse"),
    sampler = sampler,
    n_permutations = 2L,
    n_samples = 20L
  )

  sage$compute()
  result = sage$importance()

  expect_importance_dt(result, features = sage$features)
  # All importance values should be essentially zero
  checkmate::expect_numeric(result$importance, lower = -1e-10, upper = 1e-10)
})

# -----------------------------------------------------------------------------
# Sensible results
# -----------------------------------------------------------------------------

test_that("ConditionalSAGE produces sensible ranking", {
  # Use sim_dgp_independent (5 features) with Gaussian sampler
  task = sim_dgp_independent(n = 200)
  sampler = ConditionalGaussianSampler$new(task)

  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    measure = msr("regr.mse"),
    sampler = sampler,
    n_permutations = 2L,
    n_samples = 20L
  )

  sage$compute()
  result = sage$importance()
  expect_importance_dt(result, features = sage$features)

  # Important features should have higher mean importance than unimportant
  important_scores = result[grepl("^important", feature)]$importance
  unimportant_scores = result[grepl("^unimportant", feature)]$importance
  expect_gt(mean(important_scores), mean(unimportant_scores))
})

# -----------------------------------------------------------------------------
# Sampler behavior
# -----------------------------------------------------------------------------

test_that("ConditionalSAGE uses ConditionalARFSampler by default", {
  skip_if_not_installed("arf")

  # Just check default sampler class - no compute needed
  task = sim_dgp_correlated(n = 50)

  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    n_permutations = 2L,
    n_samples = 20L
  )

  checkmate::expect_r6(sage$sampler, "ConditionalARFSampler")
  expect_equal(sage$label, "Conditional SAGE")
})

test_that("ConditionalSAGE with ARF sampler computes correctly", {
  skip_on_cran() # ARF sampling is slow
  skip_if_not_installed("arf")

  # Use small task for ARF test
  task = sim_dgp_correlated(n = 50)

  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    n_permutations = 2L,
    n_samples = 20L
  )

  sage$compute()
  expect_importance_dt(sage$importance(), features = sage$features)
})

test_that("ConditionalSAGE with custom ARF sampler settings", {
  skip_on_cran() # ARF sampling is slow
  skip_if_not_installed("arf")

  task = sim_dgp_correlated(n = 50)
  custom_sampler = ConditionalARFSampler$new(task, finite_bounds = "local")

  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    measure = msr("regr.mse"),
    sampler = custom_sampler,
    n_permutations = 2L,
    n_samples = 20L
  )

  checkmate::expect_r6(sage$sampler, "ConditionalSampler")
  sage$compute()
  expect_importance_dt(sage$importance(), features = sage$features)
})

# -----------------------------------------------------------------------------
# Parameter validation
# -----------------------------------------------------------------------------

test_that("ConditionalSAGE requires predict_type='prob' for classification", {
  task = tgen("2dnormals")$generate(n = 50)

  expect_error(
    ConditionalSAGE$new(
      task = task,
      learner = lrn("classif.rpart", predict_type = "response"),
      sampler = ConditionalGaussianSampler$new(task)
    ),
    "Classification learners require probability predictions for SAGE."
  )
})

# -----------------------------------------------------------------------------
# Batching (skip_on_cran - tests implementation detail)
# -----------------------------------------------------------------------------

test_that("ConditionalSAGE batching produces consistent results", {
  skip_on_cran() # tests implementation detail, not core functionality
  skip_if_not_installed("withr")

  # Use small task with Gaussian sampler for faster batching test
  task = sim_dgp_correlated(n = 50)
  sampler = ConditionalGaussianSampler$new(task)

  # Results should be identical with or without batching
  result_batch = withr::with_seed(42, {
    sage = ConditionalSAGE$new(
      task = task,
      learner = lrn("regr.rpart"),
      sampler = ConditionalGaussianSampler$new(task),
      n_permutations = 2L,
      n_samples = 20L
    )
    sage$compute(batch_size = 1)
    sage$importance()
  })

  result_normal = withr::with_seed(42, {
    sage = ConditionalSAGE$new(
      task = task,
      learner = lrn("regr.rpart"),
      sampler = ConditionalGaussianSampler$new(task),
      n_permutations = 2L,
      n_samples = 20L
    )
    sage$compute()
    sage$importance()
  })

  expect_equal(result_batch$importance, result_normal$importance, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# n_samples parameter
# -----------------------------------------------------------------------------

test_that("ConditionalSAGE with custom n_samples", {
  # Use small task with Gaussian sampler
  task = sim_dgp_correlated(n = 50)
  sampler = ConditionalGaussianSampler$new(task)

  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    sampler = sampler,
    n_permutations = 2L,
    n_samples = 20L
  )
  expect_equal(sage$param_set$values$n_samples, 20L)

  sage$compute()
  expect_importance_dt(sage$importance(), features = sage$features)
})

# -----------------------------------------------------------------------------
# Convergence tracking (skip_on_cran - tests advanced feature)
# -----------------------------------------------------------------------------

test_that("ConditionalSAGE SE tracking in convergence_history", {
  skip_on_cran() # tests convergence tracking feature, not core SAGE

  # Use small task with Gaussian sampler
  task = sim_dgp_correlated(n = 50)
  sampler = ConditionalGaussianSampler$new(task)

  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    measure = msr("regr.mse"),
    sampler = sampler,
    n_permutations = 6L,
    n_samples = 20L
  )

  # The convergence history is tracked regardless of early stopping
  sage$compute(check_interval = 2L)

  # Check that convergence_history exists and has SE column
  expect_false(is.null(sage$convergence_history))
  expect_contains(colnames(sage$convergence_history), "se")

  # Check structure of convergence_history
  expected_cols = c("budget", "n_evals", "feature", "importance", "se")
  expect_setequal(colnames(sage$convergence_history), expected_cols)

  # SE values should be non-negative and finite
  se_values = sage$convergence_history$se
  checkmate::expect_numeric(se_values[sage$convergence_history$budget > 1], lower = 0, finite = TRUE)

  # All features should be represented in convergence history
  expect_setequal(
    unique(sage$convergence_history$feature),
    sage$features
  )
})

# -----------------------------------------------------------------------------
# Regression guard for the samples_per_row migration (+ underlying ARF fix)
# -----------------------------------------------------------------------------

test_that("ConditionalSAGE noise feature receives ~0 importance with ARFSampler", {
  skip_if_not_installed("arf")
  skip_if_not_installed("ranger")
  skip_if_not_installed("mlr3learners")

  # `sim_dgp_correlated`: y = 2 * x1 + x3 + noise; x1 correlated with x2 (r),
  # x4 is pure noise (independent, no causal effect). The expected SAGE pattern:
  #   - x1 dominant (signal, coef 2)
  #   - x3 next (signal, coef 1)
  #   - x2 some attribution via Shapley (correlated with x1 but no causal effect)
  #   - x4 ~ 0 (pure noise)
  #
  # A regression in the multi-draw conditional sampling path tends to inflate
  # all SAGE values (notably for x4) because empty-coalition predictions become
  # garbage. This test pins both the signal/noise pattern and concrete values
  # captured on the post-fix branch.
  set.seed(42L)
  task = sim_dgp_correlated(n = 200L, r = 0.7)
  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.ranger", num.trees = 100L),
    measure = msr("regr.mse"),
    sampler = ConditionalARFSampler$new(task, verbose = FALSE),
    n_permutations = 30L,
    n_samples = 30L
  )
  sage$compute()
  imp = sage$importance()

  # Correctness: noise feature x4 stays near zero, signal features dominate.
  expect_lt(abs(imp[feature == "x4", importance]), 0.3)
  expect_gt(imp[feature == "x1", importance], 0.5)
  expect_gt(imp[feature == "x1", importance], imp[feature == "x4", importance] + 1.0)

  # Regression guard. Captured on the post-fix branch at the same settings.
  # Tolerance generous enough to absorb run-to-run MC variation under different
  # R / arf / ranger versions; if it ever breaks, regenerate by running the same
  # block and pasting `imp$importance` below.
  expected_features = c("x1", "x2", "x3", "x4")
  expected_importance = c(
    1.75718693344525,
    0.359611168661264,
    0.651142439812633,
    -0.0157528350794578
  )
  expect_equal(imp$feature, expected_features)
  expect_equal(imp$importance, expected_importance, tolerance = 0.5)
})
