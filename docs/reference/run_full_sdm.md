# Run an SDM using no cross validation and previously established tune arguments

Run an SDM using no cross validation and previously established tune
arguments

## Usage

``` r
run_full_sdm(
  prep,
  tune,
  out_dir,
  return_val = "path",
  use_metric = "combo",
  force_new = FALSE,
  do_gc = FALSE,
  ...
)
```

## Arguments

- prep:

  Character or named list. If character, the path to an existing
  `prep.rds`. Otherwise, the result of a call to prep_sdm with
  return_val = "object"

- tune:

  Character or named list. If character, the path to an existing
  `tune.rds`. Otherwise, the result of a call to tune_sdm with
  return_val = "object"

- out_dir:

  FALSE or character. If FALSE the result of `run_full_sdm()` will be
  saved to a temporary folder. If character, a file 'tune.rds' will be
  created at the path defined by `out_dir`.

- return_val:

  Character: "object" or "path". Both return a named list. In the case
  of "path" the named list is simply list(full_run = out_dir). Will be
  set to "object" if `out_dir` is FALSE.

- use_metric:

  Character. Which metric to use to find the 'best' tune arguments from
  previous tuning results? Default is `combo`, the product of `auc_po`,
  `CBI_rescale` and `IMAE`. `use_metric` must be `combo` or have been
  used in the use_metrics argument to
  [`tune_sdm()`](https://acanthiza.github.io/envSDM/reference/tune_sdm.md).

- force_new:

  Logical. If outputs already exist, should they be remade?

- do_gc:

  Logical. Run `base::rm(list = ls)` and
  [`base::gc()`](https://rdrr.io/r/base/gc.html) at end of function?
  Useful when running SDMs for many, many taxa, especially if done in
  parallel.

- ...:

  Passed to
  [`tune_sdm()`](https://acanthiza.github.io/envSDM/reference/tune_sdm.md)

## Value

If `return_val` is "object" a named list. If `return_val` is "path" a
path to the saved file. If `out_dir` is a valid path, the 'full result'
(irrespective of `return_val`) is also saved to
`fs::path(out_dir, "prep.rds")`. The 'full result' is a named list with
elements:

## Examples

``` r
inst/examples/predict_sdm_ex.R
#> Error: object 'inst' not found
```
