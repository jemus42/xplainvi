# Exact SAGE estimator (estimator = "exact"): enumerates all 2^p coalitions and computes
# the Shapley decomposition in closed form, so it has no coalition-sampling error.

test_that("MarginalSAGE exact estimator works for regression and classification", {
  set.seed(3163)
  task_regr = sim_dgp_independent(n = 150)
  sage_regr = MarginalSAGE$new(
    task = task_regr,
    learner = lrn("regr.rpart"),
    estimator = "exact",
    n_samples = 20L
  )
  expect_identical(sage_regr$param_set$values$estimator, "exact")
  expect_null(sage_regr$param_set$values$n_permutations)
  expect_identical(sage_regr$converged, NA)
  sage_regr$compute()
  expect_importance_dt(sage_regr$importance(), features = sage_regr$features)
  # No coalition-sampling error, so recomputation is deterministic.
  imp1 = sage_regr$importance()$importance
  sage_regr$compute()
  expect_equal(sage_regr$importance()$importance, imp1, tolerance = 1e-12)

  budget = sage_regr$budget
  expect_equal(budget$estimator, "exact")
  expect_equal(budget$unit, "coalitions")
  expect_equal(budget$requested, 2^length(sage_regr$features))
  expect_equal(budget$used, budget$requested)
  expect_equal(budget$n_evals, budget$requested)
  expect_identical(budget$converged, NA)

  task_binary = tgen("2dnormals")$generate(n = 100)
  sage_binary = MarginalSAGE$new(
    task = task_binary,
    learner = lrn("classif.rpart", predict_type = "prob"),
    estimator = "exact",
    n_samples = 20L
  )
  sage_binary$compute()
  expect_importance_dt(sage_binary$importance(), features = sage_binary$features)

  task_multi = tgen("cassini")$generate(n = 100)
  sage_multi = MarginalSAGE$new(
    task = task_multi,
    learner = lrn("classif.rpart", predict_type = "prob"),
    estimator = "exact",
    n_samples = 20L
  )
  sage_multi$compute()
  expect_importance_dt(sage_multi$importance(), features = sage_multi$features)
  expect_length(task_multi$class_names, 3L)
})

test_that("ConditionalSAGE exact estimator works with Gaussian sampler", {
  set.seed(5417)
  task = sim_dgp_correlated(n = 120)
  sage = ConditionalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    sampler = ConditionalGaussianSampler$new(task),
    estimator = "exact",
    n_samples = 20L
  )
  sage$compute()
  expect_importance_dt(sage$importance(), features = sage$features)
})

test_that("MarginalSAGE exact with a single feature returns the total", {
  task = tgen("friedman1")$generate(n = 120)
  sage = MarginalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    features = "important4",
    estimator = "exact",
    n_samples = 20L
  )
  sage$compute()
  imp = sage$importance()
  expect_importance_dt(imp, features = "important4")
  expect_equal(nrow(imp), 1L)
})

test_that("exact estimator matches an independent brute-force Shapley computation", {
  set.seed(202)
  task = tgen("friedman1")$generate(n = 200)
  task$select(c("important1", "important2", "important4", "unimportant1"))

  sage = MarginalSAGE$new(
    task = task,
    learner = lrn("regr.rpart"),
    measure = msr("regr.mse"),
    estimator = "exact",
    n_samples = 30L
  )
  sage$compute()
  est = sage$importance()[match(sage$features, feature), importance]
  phi = brute_force_shapley(sage, task)
  expect_equal(est, unname(phi), tolerance = 1e-8)
})

test_that("permutation estimator converges to the exact estimator", {
  set.seed(808)
  task = tgen("friedman1")$generate(n = 150)
  task$select(c("important1", "important2", "important4", "unimportant1"))
  learner = lrn("regr.rpart")
  measure = msr("regr.mse")
  resampling = rsmp("holdout")$instantiate(task)

  # Share the reference subsample so both estimators target the same value function.
  set.seed(1097)
  exact = MarginalSAGE$new(task, learner, measure, resampling = resampling, estimator = "exact", n_samples = 20L)
  exact$compute()
  set.seed(1097)
  perm = MarginalSAGE$new(task, learner, measure, resampling = resampling, n_permutations = 200L, n_samples = 20L)
  perm$compute()

  cmp = merge(exact$importance(), perm$importance(), by = "feature")
  # Totals match exactly (efficiency); values match up to Monte Carlo error.
  expect_equal(sum(cmp$importance.x), sum(cmp$importance.y), tolerance = 1e-8)
  se = perm$convergence_history[budget == 200L][match(cmp$feature, feature), se]
  expect_lte(max(abs(cmp$importance.x - cmp$importance.y) - 4 * se), 0)
})

test_that("exact estimator guards its arguments", {
  task = tgen("friedman1")$generate(n = 60) # 10 features
  learner = lrn("regr.rpart")

  expect_error(MarginalSAGE$new(task, learner, estimator = "exact", max_features = 6L), "max_features")
  expect_no_error(MarginalSAGE$new(task, learner, estimator = "exact"))
  expect_error(MarginalSAGE$new(task, learner, estimator = "exact", n_permutations = 5L), "no sampling budget")
  expect_warning(MarginalSAGE$new(task, learner, estimator = "exact", early_stopping = TRUE), "ignored")
  expect_warning(MarginalSAGE$new(task, learner, max_features = 8L), "ignored")

  sage = MarginalSAGE$new(sim_dgp_independent(n = 80), learner, estimator = "exact", n_samples = 20L)
  expect_warning(sage$compute(se_threshold = 0.1), "ignored")
  expect_null(sage$convergence_history)
  expect_error(sage$plot_convergence(), "not applicable to the exact estimator")
  expect_error(suppressWarnings(sage$n_permutations <- 5L), "only valid")
})

test_that("permutation budget at or above exact cost is pointed out", {
  task = sim_dgp_independent(n = 60) # 5 features: 2^5 = 32 <= 1 + 7 * 5
  sage = MarginalSAGE$new(task, lrn("regr.rpart"), n_permutations = 7L, n_samples = 10L)
  old = xplain_opt(verbose = TRUE)
  withr::defer(xplain_opt(verbose = old$verbose))
  expect_message(sage$compute(), "estimator = \"exact\"")
})
