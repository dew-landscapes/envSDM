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
#> Error: No such file: /home/nwilloug/tmp/R/RtmpW54ve1/temp_libpath13645b38c42c8c/envSDM/examples/0.3__10__FALSE/full_run.rds
  algo <- full_run$tune_mean$algo[[1]]
#> Error: object 'full_run' not found
  model <- full_run[[paste0("tune_", algo)]]$m[[1]]
#> Error: object 'full_run' not found

  presences <- prep$testing[prep$testing$pa == 1, ]
  background <- prep$testing[prep$testing$pa == 0, ]

  evaluate_sdm(full_run$tune_rf$m[[1]]
               , p_test = presences
               , b_test = background
               )
#> Error: object 'full_run' not found
```
