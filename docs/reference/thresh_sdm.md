# Threshold a previously predicted SDM

Threshold a previously predicted SDM

## Usage

``` r
thresh_sdm(
  pred_file,
  this_taxa = NULL,
  threshold,
  thresh_file = NULL,
  terra_options = NULL,
  force_new = FALSE,
  do_gc = FALSE,
  check_tifs = TRUE
)
```

## Arguments

- pred_file:

  Character. File path of predicted sdm to threshold.

- this_taxa:

  Character. If left as default `NULL` an attempt will be made to
  extract a taxa name from `pred_file`

- threshold:

  Numeric. \> 0 and \< 1. Threshold to apply to the raster stored in the
  file at `pred_file`. Often this value will be available within the
  result of a call to
  [`tune_sdm()`](https://acanthiza.github.io/envSDM/reference/tune_sdm.md).
  e.g. `mod <- rio::import("tune.rds")` and then
  `mod$e[[1]]@thresholds$max_spec_sens`

- thresh_file:

  Character. Name to give the output threshold. If left as default
  `NULL`, `thresh_file` will be set to
  `gsub("pred", "thresh", pred_file)`

- terra_options:

  Passed to
  [`terra::terraOptions()`](https://rspatial.github.io/terra/reference/terraOptions.html).
  e.g. list(memfrac = 0.6)

- force_new:

  Logical. If output files already exist, should they be remade?

- do_gc:

  Logical. Run `base::rm(list = ls)` and
  [`base::gc()`](https://rdrr.io/r/base/gc.html) at end of function?
  Useful to keep RAM use down when running SDMs for many, many taxa,
  especially if done in parallel.

- check_tifs:

  Logical. Check if any output `.tif` files error on
  [`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html)
  and delete them if they do. Useful after a crash during pred_file.

## Value

Character path to threshold file, usually 'thresh.tif'. Output .tif and
.log, written to `out_dir`.

## Examples

``` r
  # setup -------
  source(fs::path("inst", "examples", "prep_sdm_ex.R")) # need 'data' object
#> Warning: cannot open file 'inst/examples/prep_sdm_ex.R': No such file or directory
#> Error in file(filename, "r", encoding = encoding): cannot open the connection

  # data ------
  extract_thresh <- function(tune, metric = "combo", thresh_type = "max_spec_sens") {

    tune |>
      dplyr::filter(!!rlang::ensym(metric) == max(!!rlang::ensym(metric))) |>
      dplyr::pull(!!rlang::ensym(thresh_type))

  }

  data <- data |>
    dplyr::mutate(tune_mean = purrr::map(full_run, \(x) rio::import(x, trust = TRUE)$tune_mean |> dplyr::select(algo, combo, tune_args, auc_po, max_spec_sens))
                  , threshold = purrr::map_dbl(tune_mean
                                               , extract_thresh
                                               )
                  , info = basename(dirname(out_dir))
                  )
#> Error in UseMethod("mutate"): no applicable method for 'mutate' applied to an object of class "function"


  ## thresh -------
  purrr::pwalk(list(data$pred
                    , data$threshold
                    , data$taxa
                    )
               , \(a, b, c) thresh_sdm(pred_file = a
                                       , threshold = b
                                       , this_taxa = c
                                       , thresh_file = "thresh.tif"
                                       #, force_new = TRUE
                                       )
               )
#> Error in data$pred: object of type 'closure' is not subsettable

  ## visualise-------

  tifs <- data$thresh
#> Error in data$thresh: object of type 'closure' is not subsettable

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
  terra::plot(r, cex.main = 0.7, nc = 2)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'plot': object 'r' not found
```
