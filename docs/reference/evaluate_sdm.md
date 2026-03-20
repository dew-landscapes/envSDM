# Evaluate an SDM

Returns various evaluation metrics from
[`predicts::pa_evaluate()`](https://rdrr.io/pkg/predicts/man/pa_evaluate.html)
and
[`flexsdm::sdm_eval()`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.html).

## Usage

``` r
evaluate_sdm(m, p_test, b_test, do_gc = FALSE, ...)
```

## Arguments

- m:

  SDM result within
  [`tune_sdm()`](https://acanthiza.github.io/envSDM/reference/tune_sdm.md)

- p_test:

  Presence test data generated within
  [`tune_sdm()`](https://acanthiza.github.io/envSDM/reference/tune_sdm.md)

- b_test:

  Background test data generated within
  [`tune_sdm()`](https://acanthiza.github.io/envSDM/reference/tune_sdm.md)

- do_gc:

  Logical. Run `base::rm(list = ls)` and
  [`base::gc()`](https://rdrr.io/r/base/gc.html) at end of function?
  Useful when running SDMs for many, many taxa, especially if done in
  parallel. Note, actually usees `rm(list = ls(pattern = "^[^e$]"))`.

- ...:

  Passed to both
  [`terra::predict()`](https://rspatial.github.io/terra/reference/predict.html)
  and
  [`predicts::pa_evaluate()`](https://rdrr.io/pkg/predicts/man/pa_evaluate.html)

## Value

paModelEvaluation (see
[`predicts::pa_evaluate()`](https://rdrr.io/pkg/predicts/man/pa_evaluate.html))
with extra metrics from
[`flexsdm::sdm_eval()`](https://sjevelazco.github.io/flexsdm/reference/sdm_eval.html):
AUC as auc_po_flexsdm; BOYCE as CBI; CBI_rescale (CBI is -1 to 1,
CBI_rescale is 0 to 1); and IMAE.

## Examples

``` r
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  preps <- fs::dir_ls(out_dir, regexp = "prep.rds", recurse = TRUE)

  prep <- rio::import(preps[[1]], trust = TRUE)

  full_run <- rio::import(fs::path(dirname(preps[[1]]), "combo", "full_run.rds"), trust = TRUE)
  algo <- full_run$tune_mean$algo[[1]]
  model <- full_run[[paste0("tune_", algo)]]$m[[1]]

  presences <- prep$testing[prep$testing$pa == 1, ]
  background <- prep$testing[prep$testing$pa == 0, ]

  evaluate_sdm(full_run$tune_rf$m[[1]]
               , p_test = presences
               , b_test = background
               )
#> @stats
#>   np   na prevalence   auc   cor pcor   ODP auc_po auc_po_flexsdm   CBI
#> 1 92 1000      0.084 0.991 0.738    0 0.916  0.991          0.991 0.808
#>   CBI_rescale  IMAE
#> 1       0.904 0.862
#> 
#> @thresholds
#>   max_kappa max_spec_sens no_omission equal_prevalence equal_sens_spec or10
#> 1     0.902         0.552       0.552            0.084           0.652 0.72
#> 
#> @tr_stats
#>     treshold kappa  CCR  TPR  TNR  FPR  FNR  PPP  NPP  MCR     OR
#> 1          0     0 0.08    1    0    1    0 0.08  NaN 0.92    NaN
#> 2          0  0.03 0.21    1 0.13 0.87    0  0.1    1 0.79    Inf
#> 3          0  0.04 0.27    1 0.21 0.79    0  0.1    1 0.73    Inf
#> 4        ...   ...  ...  ...  ...  ...  ...  ...  ...  ...    ...
#> 282        1  0.29 0.93 0.18    1    0 0.82 0.89 0.93 0.07 113.11
#> 283        1  0.29 0.93 0.18    1    0 0.82 0.89 0.93 0.07 113.11
#> 284        1     0 0.92    0    1    0    1  NaN 0.92 0.08    NaN
```
