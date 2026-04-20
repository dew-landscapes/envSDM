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

full_run <- rio::import(fs::path(dirname(preps[[1]]), "full_run.rds"), trust = TRUE)
algo <- full_run$tune_mean$algo[[1]]
model <- full_run[[paste0("tune_", algo)]]$m[[1]]

presences <- prep$testing[prep$testing$pa == 1, ]
background <- prep$testing[prep$testing$pa == 0, ]

evaluate_sdm(full_run$tune_rf$m[[1]]
             , p_test = presences
             , b_test = background
             )
#> @stats
#>   np  na prevalence   auc   cor pcor   ODP auc_po auc_po_flexsdm   CBI
#> 1 97 184      0.345 0.973 0.888    0 0.655  0.973          0.973 0.911
#>   CBI_rescale  IMAE
#> 1       0.955 0.901
#> 
#> @thresholds
#>   max_kappa max_spec_sens no_omission equal_prevalence equal_sens_spec  or10
#> 1     0.742         0.742       0.032            0.346           0.414 0.822
#> 
#> @tr_stats
#>     treshold kappa  CCR  TPR  TNR  FPR  FNR  PPP  NPP  MCR     OR
#> 1          0     0 0.35    1    0    1    0 0.35  NaN 0.65    NaN
#> 2          0  0.12 0.45    1 0.16 0.84    0 0.39    1 0.55    Inf
#> 3          0  0.18  0.5    1 0.24 0.76    0 0.41    1  0.5    Inf
#> 4        ...   ...  ...  ...  ...  ...  ...  ...  ...  ...    ...
#> 101        1  0.76  0.9 0.71 0.99 0.01 0.29 0.99 0.87  0.1 450.96
#> 102        1  0.76  0.9 0.71 0.99 0.01 0.29 0.99 0.87  0.1 450.96
#> 103        1     0 0.65    0    1    0    1  NaN 0.65 0.35    NaN
```
