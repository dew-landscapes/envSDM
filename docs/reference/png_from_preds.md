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
out_dir <- file.path(system.file(package = "envSDM"), "examples")

# setup -------
data <- readRDS(fs::path(out_dir, "data.rds"))

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
#> [[1]]
#> /projects/dev/nige/packages/envSDM/inst/examples/0__10__TRUE/full_run.rds
#> 
#> [[2]]
#> /projects/dev/nige/packages/envSDM/inst/examples/0.3__10__TRUE/full_run.rds
#> 
#> [[3]]
#> /projects/dev/nige/packages/envSDM/inst/examples/0__30__TRUE/full_run.rds
#> 
#> [[4]]
#> /projects/dev/nige/packages/envSDM/inst/examples/0.3__30__TRUE/full_run.rds
#> 
#> [[5]]
#> /projects/dev/nige/packages/envSDM/inst/examples/0__10__FALSE/full_run.rds
#> 
#> [[6]]
#> /projects/dev/nige/packages/envSDM/inst/examples/0.3__10__FALSE/full_run.rds
#> 
#> [[7]]
#> /projects/dev/nige/packages/envSDM/inst/examples/0__30__FALSE/full_run.rds
#> 
#> [[8]]
#> /projects/dev/nige/packages/envSDM/inst/examples/0.3__30__FALSE/full_run.rds
#> 


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
#> Error in pmap(.l, .f, ..., .progress = .progress): ℹ In index: 1.
#> Caused by error in `.f()`:
#> ! object 'preds' not found

## visualise-------
tifs <- fs::path(data$out_dir[file.exists(fs::path(data$out_dir, "pred.tif"))], "pred.tif")

names <- paste0("hold_prop "
                , data$hold_prop
                , "; stretch "
                , data$stretch
                , "; new_bg "
                , data$new_bg_test
                )

r <- terra::rast(tifs)
names(r) <- names
terra::panel(r, cex.main = 0.6, nc = 2)
```
