# Generate .png (map) files for each prediction

Finds any .tif files in `pred_dir` and writes them to .png files.
Includes the retrieval and addition to the map of: various SDM metrics;
and the original presence points.

## Usage

``` r
png_from_preds(
  pred_dir,
  full_run_dir = NULL,
  trim = TRUE,
  force_new = FALSE,
  do_gc = TRUE,
  ...
)
```

## Arguments

- trim:

  Logical. Trim NA values from outside (using
  [`terra::trim()`](https://rspatial.github.io/terra/reference/trim.html))

- force_new:

  Logical. If .png file already exists, recreate it?

- do_gc:

  Logical. Run `base::rm(list = ls)` and
  [`base::gc()`](https://rdrr.io/r/base/gc.html) at end of function?
  Useful when running SDMs for many, many taxa, especially if done in
  parallel.

- ...:

  Passed to [`fs::dir_ls()`](https://fs.r-lib.org/reference/dir_ls.html)

- prep:

  Character or named list. If character, the path to an existing
  `prep.rds`. Otherwise, the result of a call to
  [`prep_sdm()`](https://acanthiza.github.io/envSDM/reference/prep_sdm.md)
  with return_val = "object".

- full_run:

  Character or named list. If character, the path to an existing
  `full_run.rds`. Otherwise, the result of a call to `run_full_sdm`()
  with return_val = "object".

- out_dir:

  Character. Name of directory into which `.pngs`s will be saved. Will
  be created if it does not exist.

## Value

`invisible(NULL)`. Writes .png files with the same file name as any .tif
files

## Examples

``` r
  # setup -------
  source(fs::path("inst", "examples", "prep_sdm_ex.R")) # need 'data' object
#> Warning: cannot open file 'inst/examples/prep_sdm_ex.R': No such file or directory
#> Error in file(filename, "r", encoding = encoding): cannot open the connection

  # Best combo--------
  ## run full SDM --------
  purrr::pmap(list(data$prep
                    , data$tune
                    , data$out_dir
                    )
               , \(a, b, c) run_full_sdm(prep = a
                                         , tune = b
                                         , out_dir = c
                                         , use_metric = "combo"

                                         # passed to tune_sdm via dots
                                         , metrics_df = envSDM::sdm_metrics
                                         #, force_new = FALSE
                                         )
               )
#> Error in data$prep: object of type 'closure' is not subsettable


  ## predict -------
  purrr::pwalk(list(data$prep
                    , data$out_dir
                    )
               , \(a, b) predict_sdm(prep = a
                                     , full_run = fs::path(b, "full_run.rds")
                                     , out_dir = b
                                     , predictors = preds
                                     , check_tifs = TRUE
                                     #, force_new = FALSE
                                     )
               )
#> Error in data$prep: object of type 'closure' is not subsettable

  ## visualise-------
  tifs <- fs::path(data$out_dir[file.exists(fs::path(data$out_dir, "pred.tif"))], "pred.tif")
#> Error in data$out_dir: object of type 'closure' is not subsettable

  names <- paste0("hold_prop "
                  , data$hold_prop
                  , "; stretch "
                  , data$stretch
                  , "; new_bg "
                  , data$new_bg_test
                  )
#> Error in data$hold_prop: object of type 'closure' is not subsettable

  r <- terra::rast(tifs)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rast': object 'tifs' not found
  names(r) <- names
#> Error: object 'r' not found
  terra::panel(r, cex.main = 0.6, nc = 2)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'panel': object 'r' not found
```
