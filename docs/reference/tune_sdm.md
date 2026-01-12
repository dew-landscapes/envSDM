# Tune, and evaluate, species distribution models

Tune, and evaluate, species distribution models

## Usage

``` r
tune_sdm(
  prep,
  out_dir = FALSE,
  return_val = "path",
  algo = c("all", "maxnet", "envelope", "rf"),
  max_corr = list(maxnet = 0.7, envelope = 0.9, rf = 0.99),
  fc = "auto_feature",
  limit_p = FALSE,
  rm = seq(1, 6, 0.5),
  trees = c(999),
  mtry = TRUE,
  limit_spat_mtry = 4,
  nodesize = c(1, 2),
  keep_model = FALSE,
  best_run = FALSE,
  metrics_df = envSDM::sdm_metrics,
  use_metrics = c("auc_po", "CBI_rescale", "IMAE"),
  do_gc = FALSE,
  force_new = FALSE,
  ...
)
```

## Arguments

- prep:

  Character or named list. If character, the path to an existing
  `prep.rds`. Otherwise, the result of a call to prep_sdm with
  return_val = "object"

- out_dir:

  FALSE or character. If FALSE the result of tune_sdm will be saved to a
  temporary folder. If character, a file 'tune.rds' will be created at
  the path defined by out_dir.

- return_val:

  Character: "object" or "path". Both return a named list. In the case
  of "path" the named list is simply list(tune = out_dir). Will be set
  to "object" if `out_dir` is FALSE.

- algo:

  Character. Name of algorithm to use.

- max_corr:

  Named list. Names of list elements must match algorithms being used.
  For each pair of predictor variables correlated at or above `max_corr`
  one will be dropped using
  [`caret::findCorrelation()`](https://rdrr.io/pkg/caret/man/findCorrelation.html).

- fc:

  Character. Used to generate levels of `classes` argument to
  [`maxnet::maxnet()`](https://rdrr.io/pkg/maxnet/man/maxnet.html) that
  are tuned.

- limit_p:

  `TRUE`, `FALSE` or number of predictor variables above which to limit
  the use of `p` in the classes argument used in
  [`maxnet::maxnet()`](https://rdrr.io/pkg/maxnet/man/maxnet.html).
  Useful with many predictor variables when it becomes unwieldy to
  generate interactions for all predictors.

- rm:

  Numeric. Used to generate levels of `regmult` argument to
  [`maxnet::maxnet()`](https://rdrr.io/pkg/maxnet/man/maxnet.html) that
  are tuned.

- trees:

  Used to generate the levels of `ntree` argument to
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
  that are tuned. `TRUE` (tune with default `trees`), `FALSE` (don't
  tune `trees`) or numeric (the `trees` values to tune with).

- mtry:

  Used to generate the levels of `mtry` argument to
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
  that are tuned. `TRUE` (tune with sensible guesses for `mtry`),
  `FALSE` (only use default
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
  `mtry`) or numeric (the `mtry` values to tune with).

- limit_spat_mtry:

  Numeric. If `mtry` is `TRUE` and if using spatial cross validation,
  the values of `mtry` to tune will be limited to less than or equal to
  `limit_spat_mtry`.

- nodesize:

  Used to generate the levels of `nodesize` argument to
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
  that are tuned. `TRUE` (tune with default `nodesize`), `FALSE` (only
  use default
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
  `nodesize`) or numeric (the `nodesize` values to tune with).

- keep_model:

  Logical. If `TRUE` the model results will be appended as a list column
  in the returned tibble (as column `m`)

- best_run:

  Logical. If `TRUE` this alters the behaviour of the `tune_sdm()` by,
  well, not tuning. :). Sets all blocks to the same value so no
  cross-validation.

- metrics_df:

  Dataframe. Defines which metrics to use when deciding on 'good' SDMs.

- use_metrics:

  Character. Vector of values in metrics_df\$metric to use when finding
  the 'best' model.

- do_gc:

  Logical. Run `base::rm(list = ls)` and
  [`base::gc()`](https://rdrr.io/r/base/gc.html) at end of function?
  Useful when running SDMs for many, many taxa, especially if done in
  parallel.

- force_new:

  Logical. If outputs already exist, should they be remade?

- ...:

  Passed to
  [`evaluate_sdm()`](https://acanthiza.github.io/envSDM/reference/evaluate_sdm.md).
  e.g. thresholds for use in
  [`predicts::pa_evaluate()`](https://rdrr.io/pkg/predicts/man/pa_evaluate.html)
  (as `tr` argument, although if used, the values of the `thresholds`
  element of the `pa_ModelEvaluation` object returned by
  [`predicts::pa_evaluate()`](https://rdrr.io/pkg/predicts/man/pa_evaluate.html)
  will be limited to the values in `tr`).

## Value

If `return_val` is "object" a named list. If `return_val` is "path" a
path to the saved file. If `out_dir` is a valid path, the 'full result'
(irrespective of `return_val`) is also saved to
`fs::path(out_dir, "prep.rds")`. The 'full result' is a named list with
elements:

## Examples

``` r
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  data <- fs::path(system.file(package = "envSDM"), "examples") |>
    fs::dir_ls(regexp = "prep\\.rds$"
               , recurse = TRUE
               ) |>
    tibble::enframe(name = NULL, value = "prep") |>
    dplyr::mutate(taxa = gsub("\\.rds", "", basename(dirname(prep)))
                  , out_dir = fs::path(out_dir, taxa)
                  )

  purrr::map(data$out_dir
              , \(x) tune_sdm(prep = fs::path(x, "prep.rds")
                              , out_dir = x
                              , fc = "lq"
                              , rm = c(2, 3)
                              , trees = 500
                              , mtry = c(1:3)
                              , nodesize = c(1, 2, 3)
                              , limit_p = 3
                              , use_metrics = c("auc_po", "CBI_rescale", "IMAE", "or10")
                              #, force_new = TRUE
                              )
              )
#> [[1]]
#> /home/nwilloug/temp/RtmpZmDzsp/temp_libpath19cb0544a4dd6f/envSDM/examples/chg/tune.rds
#> 
#> [[2]]
#> /home/nwilloug/temp/RtmpZmDzsp/temp_libpath19cb0544a4dd6f/envSDM/examples/mjs/tune.rds
#> 
#> [[3]]
#> /home/nwilloug/temp/RtmpZmDzsp/temp_libpath19cb0544a4dd6f/envSDM/examples/wjb/tune.rds
#> 

  # which tune args were best for each taxa using 'combo'?
  data %>%
    dplyr::mutate(tune = fs::path(out_dir, "tune.rds")
                  , tune = purrr::map(tune, rio::import, trust = TRUE)
                  , tune_mean = purrr::map(tune, "tune_mean")
                  ) %>%
    tidyr::unnest(cols = c(tune_mean)) %>%
    dplyr::filter(best) %>% # used 'combo' to determine 'best' as default in tune_sdm
    dplyr::select(taxa, algo, tune_args, combo, auc_po, IMAE, CBI, max_spec_sens)
#> # A tibble: 0 × 8
#> # ℹ 8 variables: taxa <chr>, algo <chr>, tune_args <chr>, combo <dbl>,
#> #   auc_po <dbl>, IMAE <dbl>, CBI <dbl>, max_spec_sens <dbl>

  # or best tune args choosing on just auc_po?
  data %>%
    dplyr::mutate(tune = fs::path(out_dir, "tune.rds")
                  , tune = purrr::map(tune, rio::import, trust = TRUE)
                  , all = purrr::map(tune, "tune_mean")
                  ) %>%
    tidyr::unnest(cols = c(all)) %>%
    dplyr::group_by(taxa) %>%
    dplyr::filter(auc_po == max(auc_po)) %>%
    dplyr::ungroup() %>%
    dplyr::select(taxa, algo, tune_args, auc_po, IMAE, CBI, max_spec_sens)
#> # A tibble: 2 × 7
#>   taxa  algo   tune_args             auc_po  IMAE   CBI max_spec_sens
#>   <chr> <chr>  <chr>                  <dbl> <dbl> <dbl>         <dbl>
#> 1 chg   rf     tr: 500. mt: 2. ns: 2  0.772 0.812 0.394         0.276
#> 2 mjs   maxnet fc: lq. rm: 3          0.832 0.840 0.629         0.258
```
