# Shapley Additive Global Importance (SAGE) Base Class

Base class for SAGE (Shapley Additive Global Importance) feature
importance based on Shapley values with marginalization. This is an
abstract class - use
[MarginalSAGE](https://mlr-org.github.io/xplainfi/reference/MarginalSAGE.md)
or
[ConditionalSAGE](https://mlr-org.github.io/xplainfi/reference/ConditionalSAGE.md).

## Details

SAGE uses Shapley values to fairly distribute the total prediction
performance among all features. Unlike perturbation-based methods, SAGE
marginalizes features by integrating over their distribution. This is
approximated by averaging predictions over a reference dataset.

SAGE values are reductions in the measure's score relative to the empty
coalition, `score(empty) - score(S)`, so that positive values mean the
feature improves performance. For measures that are maximized
(`measure$minimize = FALSE`, e.g. `classif.acc`) the scores are negated
internally, so the sign convention is the same for all measures.

**Standard Error Calculation**: The standard errors (SE) reported in
`$convergence_history` reflect the uncertainty in Shapley value
estimation across different random permutations within a single
resampling iteration. These SEs quantify the Monte Carlo sampling error
for a fixed trained model and are only valid for inference about the
importance of features for that specific model. They do not capture
broader uncertainty from model variability across different train/test
splits or resampling iterations.

## References

Covert I, Lundberg S, Lee S (2020). “Understanding Global Feature
Contributions With Additive Importance Measures.” In *Advances in Neural
Information Processing Systems*, volume 33, 17212–17223.
<https://proceedings.neurips.cc/paper/2020/hash/c7bf0b7c1a86d5eb3be2c722cf2cf746-Abstract.html>.

## See also

[MarginalSAGE](https://mlr-org.github.io/xplainfi/reference/MarginalSAGE.md)
[ConditionalSAGE](https://mlr-org.github.io/xplainfi/reference/ConditionalSAGE.md)

## Super class

[`FeatureImportanceMethod`](https://mlr-org.github.io/xplainfi/reference/FeatureImportanceMethod.md)
-\> `SAGE`

## Public fields

- `n_permutations`:

  (`integer(1)`) Number of permutations to sample.

- `convergence_history`:

  ([`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html))
  History of SAGE values during computation.

- `converged`:

  (`logical(1)`) Whether convergence was detected.

- `n_permutations_used`:

  (`integer(1)`) Actual number of permutations used.

## Methods

### Public methods

- [`SAGE$new()`](#method-SAGE-initialize)

- [`SAGE$compute()`](#method-SAGE-compute)

- [`SAGE$plot_convergence()`](#method-SAGE-plot_convergence)

- [`SAGE$clone()`](#method-SAGE-clone)

Inherited methods

- [`FeatureImportanceMethod$importance()`](https://mlr-org.github.io/xplainfi/reference/FeatureImportanceMethod.html#method-importance)
- [`FeatureImportanceMethod$obs_loss()`](https://mlr-org.github.io/xplainfi/reference/FeatureImportanceMethod.html#method-obs_loss)
- [`FeatureImportanceMethod$print()`](https://mlr-org.github.io/xplainfi/reference/FeatureImportanceMethod.html#method-print)
- [`FeatureImportanceMethod$reset()`](https://mlr-org.github.io/xplainfi/reference/FeatureImportanceMethod.html#method-reset)
- [`FeatureImportanceMethod$scores()`](https://mlr-org.github.io/xplainfi/reference/FeatureImportanceMethod.html#method-scores)

------------------------------------------------------------------------

### `SAGE$new()`

Creates a new instance of the SAGE class.

#### Usage

    SAGE$new(
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
    )

#### Arguments

- `task, learner, measure, resampling, features`:

  Passed to FeatureImportanceMethod.

- `n_permutations`:

  (`integer(1)`: `10L`) Number of permutations to sample for SAGE value
  estimation. The total number of evaluated coalitions is
  `1 (empty) + n_permutations * n_features`.

- `batch_size`:

  (`integer(1)`: `5000L`) Maximum number of observations to process in a
  single prediction call.

- `n_samples`:

  (`integer(1)`: `100L`) Number of samples to use for marginalizing
  out-of-coalition features. For
  [MarginalSAGE](https://mlr-org.github.io/xplainfi/reference/MarginalSAGE.md),
  this is the number of marginal data samples ("background data" in
  other implementations). For
  [ConditionalSAGE](https://mlr-org.github.io/xplainfi/reference/ConditionalSAGE.md),
  this is the number of conditional samples per test instance retrieved
  from `sampler`.

- `early_stopping`:

  (`logical(1)`: `TRUE`) Whether to enable early stopping based on
  convergence detection.

- `se_threshold`:

  (`numeric(1)`: `0.01`) Convergence threshold for relative standard
  error. Convergence is detected when the maximum relative SE across all
  features falls below this threshold. Relative SE is calculated as SE
  divided by the range of importance values (max - min), making it
  scale-invariant across different loss metrics. Default of `0.01` means
  convergence when relative SE is below 1% of the importance range.

- `min_permutations`:

  (`integer(1)`: `10L`) Minimum permutations before checking for
  convergence. Convergence is judged based on the standard errors of the
  estimated SAGE values, which requires a sufficiently large number of
  samples (i.e., evaluated coalitions).

- `check_interval`:

  (`integer(1)`: `1L`) Check convergence every N permutations.

------------------------------------------------------------------------

### `SAGE$compute()`

Compute SAGE values.

#### Usage

    SAGE$compute(
      store_backends = TRUE,
      batch_size = NULL,
      early_stopping = NULL,
      se_threshold = NULL,
      min_permutations = NULL,
      check_interval = NULL
    )

#### Arguments

- `store_backends`:

  (`logical(1)`) Whether to store data backends.

- `batch_size`:

  (`integer(1)`: `5000L`) Maximum number of observations to process in a
  single prediction call.

- `early_stopping`:

  (`logical(1)`: `TRUE`) Whether to check for convergence and stop
  early.

- `se_threshold`:

  (`numeric(1)`: `0.01`) Convergence threshold for relative standard
  error. SE is normalized by the range of importance values (max - min)
  to make convergence detection scale-invariant. Default `0.01` means
  convergence when relative SE \< 1%.

- `min_permutations`:

  (`integer(1)`: `10L`) Minimum permutations before checking
  convergence.

- `check_interval`:

  (`integer(1)`: `1L`) Check convergence every N permutations.

------------------------------------------------------------------------

### `SAGE$plot_convergence()`

Plot convergence history of SAGE values.

#### Usage

    SAGE$plot_convergence(features = NULL)

#### Arguments

- `features`:

  (`character` \| `NULL`) Features to plot. If NULL, plots all features.

#### Returns

A [ggplot2](https://ggplot2.tidyverse.org/reference/ggplot.html) object

------------------------------------------------------------------------

### `SAGE$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SAGE$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
