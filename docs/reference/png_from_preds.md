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
max_cores <- nrow(data)
use_cores <- min(max_cores, parallel::detectCores() - 1)

future::plan(future::multisession(workers = use_cores))
#> Error in future::plan(future::multisession(workers = use_cores)): object 'use_cores' not found

furrr::future_pwalk(list(data$prep
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


## predict -------
furrr::future_pwalk(list(data$prep
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
#> Error in (function (.l, .f, ..., .progress = FALSE) {    pmap_("list", .l, .f, ..., .progress = .progress)})(.l = list(structure(c("/projects/dev/nige/packages/envSDM/inst/examples/chg__0__10__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0.3__10__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0__30__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0.3__30__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0__10__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0.3__10__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0__30__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0.3__30__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0__10__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0.3__10__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0__30__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0.3__30__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0__10__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0.3__10__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0__30__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0.3__30__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0__10__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0.3__10__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0__30__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0.3__30__TRUE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0__10__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0.3__10__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0__30__FALSE/prep.rds", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0.3__30__FALSE/prep.rds"), class = c("fs_path", "character")), structure(c("/projects/dev/nige/packages/envSDM/inst/examples/chg__0__10__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0.3__10__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0__30__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0.3__30__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0__10__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0.3__10__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0__30__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/chg__0.3__30__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0__10__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0.3__10__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0__30__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0.3__30__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0__10__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0.3__10__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0__30__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/mjs__0.3__30__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0__10__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0.3__10__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0__30__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0.3__30__TRUE", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0__10__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0.3__10__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0__30__FALSE", "/projects/dev/nige/packages/envSDM/inst/examples/wjb__0.3__30__FALSE"), class = c("fs_path", "character"))), .f = function (...) {    NULL    NULL    ...furrr_out <- ...furrr_fn(...)    NULL}): ℹ In index: 1.
#> Caused by error in `...furrr_fn()`:
#> ! object 'preds' not found

future::plan(future::sequential())

## visualise-------
# just use one taxa
vis_data <- data |>
  dplyr::filter(taxa == "chg")

tifs <- fs::path(vis_data$out_dir[file.exists(fs::path(vis_data$out_dir, "pred.tif"))], "pred.tif")

names <- paste0("hold_prop "
                , vis_data$hold_prop
                , "; stretch "
                , vis_data$stretch
                , "; spatial_folds "
                , vis_data$spatial_folds
                )

r <- terra::rast(tifs)
names(r) <- names
terra::panel(r, cex.main = 0.6, nc = 2)
```
