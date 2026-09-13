test_that("sage_batch_predict works for regression without batching", {
  skip_if_not_installed("rpart")

  set.seed(123)
  task = tgen("friedman1")$generate(n = 50)
  learner = lrn("regr.rpart")
  learner$train(task)

  test_data = task$data()
  predictions = sage_batch_predict(learner, test_data, task, batch_size = NULL, task_type = "regr")

  expect_type(predictions, "double")
  expect_length(predictions, 50)
  expect_true(all(is.finite(predictions)))
})

test_that("sage_batch_predict works for regression with batching", {
  skip_if_not_installed("rpart")

  set.seed(123)
  task = tgen("friedman1")$generate(n = 50)
  learner = lrn("regr.rpart")
  learner$train(task)

  test_data = task$data()

  # Predict without batching
  pred_no_batch = sage_batch_predict(
    learner,
    test_data,
    task,
    batch_size = NULL,
    task_type = "regr"
  )

  # Predict with batching
  pred_with_batch = sage_batch_predict(
    learner,
    test_data,
    task,
    batch_size = 10,
    task_type = "regr"
  )

  # Should be identical
  expect_equal(pred_no_batch, pred_with_batch)
})

test_that("sage_batch_predict works for classification without batching", {
  skip_if_not_installed("rpart")

  set.seed(123)
  task = tgen("2dnormals")$generate(n = 50)
  learner = lrn("classif.rpart", predict_type = "prob")
  learner$train(task)

  test_data = task$data()
  predictions = sage_batch_predict(
    learner,
    test_data,
    task,
    batch_size = NULL,
    task_type = "classif"
  )

  expect_true(is.matrix(predictions))
  expect_equal(nrow(predictions), 50)
  expect_equal(ncol(predictions), 2) # Binary classification
  expect_true(all(is.finite(predictions)))
  # Probabilities should sum to 1
  expect_true(all(abs(rowSums(predictions) - 1) < 1e-10))
})

test_that("sage_batch_predict works for classification with batching", {
  skip_if_not_installed("rpart")

  set.seed(123)
  task = tgen("2dnormals")$generate(n = 50)
  learner = lrn("classif.rpart", predict_type = "prob")
  learner$train(task)

  test_data = task$data()

  # Predict without batching
  pred_no_batch = sage_batch_predict(
    learner,
    test_data,
    task,
    batch_size = NULL,
    task_type = "classif"
  )

  # Predict with batching
  pred_with_batch = sage_batch_predict(
    learner,
    test_data,
    task,
    batch_size = 10,
    task_type = "classif"
  )

  # Should be identical
  expect_equal(pred_no_batch, pred_with_batch)
})

test_that("sage_aggregate_predictions works for regression", {
  # Create sample data with multiple samples per coalition/instance
  combined_data = data.table::data.table(
    .coalition_id = c(1, 1, 1, 2, 2, 2),
    .test_instance_id = c(1, 1, 1, 1, 1, 1),
    feature1 = c(1, 2, 3, 4, 5, 6)
  )

  predictions = c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

  result = sage_aggregate_predictions(
    combined_data,
    predictions,
    task_type = "regr",
    class_names = NULL
  )

  expect_equal(nrow(result), 2) # 2 coalitions
  expect_true(all(c(".coalition_id", ".test_instance_id", "avg_pred") %in% names(result)))
  expect_equal(result[.coalition_id == 1]$avg_pred, 2.0) # mean(1,2,3)
  expect_equal(result[.coalition_id == 2]$avg_pred, 5.0) # mean(4,5,6)
})

test_that("sage_aggregate_predictions works for classification", {
  # Create sample data with multiple samples per coalition/instance
  combined_data = data.table::data.table(
    .coalition_id = c(1, 1, 1, 2, 2, 2),
    .test_instance_id = c(1, 1, 1, 1, 1, 1),
    feature1 = c(1, 2, 3, 4, 5, 6)
  )

  # Binary classification probabilities
  predictions = matrix(
    c(
      0.9,
      0.1,
      0.8,
      0.2,
      0.7,
      0.3,
      0.6,
      0.4,
      0.5,
      0.5,
      0.4,
      0.6
    ),
    ncol = 2,
    byrow = TRUE
  )

  class_names = c("A", "B")

  result = sage_aggregate_predictions(
    combined_data,
    predictions,
    task_type = "classif",
    class_names = class_names
  )

  expect_equal(nrow(result), 2) # 2 coalitions
  expect_true(all(c(".coalition_id", ".test_instance_id", "A", "B") %in% names(result)))

  # Check averaged probabilities for coalition 1
  expect_equal(result[.coalition_id == 1]$A, mean(c(0.9, 0.8, 0.7)))
  expect_equal(result[.coalition_id == 1]$B, mean(c(0.1, 0.2, 0.3)))

  # Check averaged probabilities for coalition 2
  expect_equal(result[.coalition_id == 2]$A, mean(c(0.6, 0.5, 0.4)))
  expect_equal(result[.coalition_id == 2]$B, mean(c(0.4, 0.5, 0.6)))
})

test_that("sage_aggregate_predictions handles multiple test instances", {
  # Create data with 2 coalitions and 2 test instances
  combined_data = data.table::data.table(
    .coalition_id = c(1, 1, 1, 1, 2, 2, 2, 2),
    .test_instance_id = c(1, 1, 2, 2, 1, 1, 2, 2),
    feature1 = 1:8
  )

  predictions = 1:8

  result = sage_aggregate_predictions(
    combined_data,
    predictions,
    task_type = "regr",
    class_names = NULL
  )

  expect_equal(nrow(result), 4) # 2 coalitions × 2 test instances
  expect_equal(result[.coalition_id == 1 & .test_instance_id == 1]$avg_pred, mean(c(1, 2)))
  expect_equal(result[.coalition_id == 1 & .test_instance_id == 2]$avg_pred, mean(c(3, 4)))
  expect_equal(result[.coalition_id == 2 & .test_instance_id == 1]$avg_pred, mean(c(5, 6)))
  expect_equal(result[.coalition_id == 2 & .test_instance_id == 2]$avg_pred, mean(c(7, 8)))
})

test_that("sage_growing_coalitions builds row-major growing prefixes", {
  perms = list(c("a", "b"), c("b", "a"))
  co = sage_growing_coalitions(perms)
  expect_length(co, 4L)
  expect_identical(co[[1]], "a")
  expect_identical(co[[2]], c("a", "b"))
  expect_identical(co[[3]], "b")
  expect_identical(co[[4]], c("b", "a"))
})

test_that("sage_marginal_contributions accumulates closed-form, offset-aware", {
  feats = c("a", "b")
  perms = list(c("a", "b"), c("b", "a"))
  baseline = 10
  # losses laid out row-major over (perm, step), no offset:
  # perm1: {a}=8, {a,b}=5 ; perm2: {b}=9, {b,a}=4
  losses = c(8, 5, 9, 4)
  r = sage_marginal_contributions(perms, losses, baseline, feats, offset = 0L)
  # perm1: a: 10-8=2 ; b: 8-5=3
  # perm2: b: 10-9=1 ; a: 9-4=5
  expect_equal(r$sv[["a"]], 2 + 5)
  expect_equal(r$sv[["b"]], 3 + 1)
  expect_equal(r$sv_sq[["a"]], 2^2 + 5^2)
  expect_equal(r$sv_sq[["b"]], 3^2 + 1^2)
  expect_setequal(names(r$sv), feats)

  # offset = 1L: an extra leading slot (e.g. empty coalition) shifts indices by 1
  losses_off = c(99, 8, 5, 9, 4)
  r2 = sage_marginal_contributions(perms, losses_off, baseline, feats, offset = 1L)
  expect_equal(r2$sv, r$sv)
  expect_equal(r2$sv_sq, r$sv_sq)
})

test_that("sage_convergence_ratio normalizes by spread and is NA without SEs", {
  expect_equal(sage_convergence_ratio(c(a = 0, b = 2, c = 1), c(0.1, 0.4, 0.2)), 0.2)
  # degenerate spread falls back to the absolute SE
  expect_equal(sage_convergence_ratio(c(1, 1), c(0.3, 0.1)), 0.3)
  expect_identical(sage_convergence_ratio(c(0, 1), c(NA_real_, 0.1)), NA_real_)
})
