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

  source(fs::path(out_dir, "tune_sdm_ex.R")) # make sure following prep file exists

  prep <- rio::import(fs::path(data$out_dir[[1]], "prep.rds"), trust = TRUE)

  model <- tune_sdm(prep = prep
                    , out_dir = FALSE
                    , return_val = "object"
                    , algo = "rf"
                    , trees = 500
                    , mtry = 2
                    , nodesize = 1
                    , keep_model = TRUE
                    )
#> tuning chg with algorithms: rf
#> out_dir is /home/nwilloug/temp/RtmpJkraMl/file19cbce2552cfb7
#> rf tune
#> tune rf ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% |  ETA:  0s

  presences <- prep$testing$testing[[1]][prep$testing$testing[[1]]$pa == 1, ]
  background <- prep$testing$testing[[1]][prep$testing$testing[[1]]$pa == 0, ]

  evaluate_sdm(model$tune_rf$m[[1]]
               , p_test = presences
               , b_test = background
               )
#> @stats
#>   np   na prevalence  auc   cor pcor   ODP auc_po auc_po_flexsdm   CBI
#> 1 92 5069      0.018 0.98 0.468    0 0.982   0.98           0.98 0.882
#>   CBI_rescale  IMAE
#> 1       0.941 0.862
#> 
#> @thresholds
#>   max_kappa max_spec_sens no_omission equal_prevalence equal_sens_spec  or10
#> 1     0.938         0.634       0.236            0.018            0.41 0.644
#> 
#> @tr_stats
#>     treshold kappa  CCR  TPR  TNR  FPR  FNR  PPP  NPP  MCR     OR
#> 1          0     0 0.02    1    0    1    0 0.02  NaN 0.98    NaN
#> 2          0     0 0.13    1 0.12 0.88    0 0.02    1 0.87    Inf
#> 3          0  0.01 0.19    1 0.17 0.83    0 0.02    1 0.81    Inf
#> 4        ...   ...  ...  ...  ...  ...  ...  ...  ...  ...    ...
#> 353        1  0.36 0.99 0.23    1    0 0.77 0.88 0.99 0.01 499.46
#> 354        1  0.36 0.99 0.23    1    0 0.77 0.88 0.99 0.01 499.46
#> 355        1     0 0.98    0    1    0    1  NaN 0.98 0.02    NaN
```
