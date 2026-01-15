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
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  # data ------
  extract_thresh <- function(tune, metric, thresh_type = "max_spec_sens") {

    tune |>
      dplyr::filter(!!rlang::ensym(metric) == max(!!rlang::ensym(metric))) |>
      dplyr::pull(!!rlang::ensym(thresh_type))

  }

  data <- fs::path(system.file(package = "envSDM"), "examples") |>
    fs::dir_ls(regexp = "pred\\.tif"
               , recurse = TRUE
               ) |>
    tibble::enframe(name = NULL, value = "pred") |>
    dplyr::mutate(out_dir = dirname(pred)
                  , taxa = basename(dirname(out_dir))
                  , metric = basename(out_dir)
                  , tune = fs::dir_ls(out_dir, regexp = "full_run.rds")
                  , tune_mean = purrr::map(tune, \(x) rio::import(x, trust = TRUE)$tune_mean |> dplyr::select(algo, tune_args, auc_po, combo, max_spec_sens))
                  , thresh = purrr::map2_dbl(tune_mean
                                             , metric
                                             , extract_thresh
                                             )

                  )

  ## thresh -------
  purrr::pwalk(list(data$pred
                    , data$thresh
                    , data$taxa
                    )
               , \(a, b, c) thresh_sdm(pred_file = a
                                       , threshold = b
                                       , this_taxa = c
                                       , thresh_file = "thresh.tif"
                                       #, force_new = TRUE
                                       )
               )
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0/auc_po/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0/combo/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3/auc_po/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3/combo/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3__1/auc_po/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3__1/combo/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3__5/auc_po/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3__5/combo/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0__1/auc_po/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0__1/combo/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0__5/auc_po/thresh.tif already exists
#> threshold file: /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0__5/combo/thresh.tif already exists

  ## visualise-------
  ### threshold -------
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "thresh.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )

















```
