# Predict from SDM

The resulting `pred.tif` is masked to the boundary provided to the
`pred_limit` argument of prep_sdm; or generated in prep_sdm from the
`pred_limit`, `limit_buffer` and `pred_clip` arguments.

## Usage

``` r
predict_sdm(
  prep,
  full_run,
  out_dir,
  file_name = "pred.tif",
  use_env_naming = FALSE,
  predictors = NULL,
  is_env_pred = FALSE,
  terra_options = NULL,
  doClamp = TRUE,
  force_new = FALSE,
  do_gc = FALSE,
  check_tifs = FALSE,
  handle_errors = TRUE,
  ...
)
```

## Arguments

- prep:

  Character or named list. If character, the path to an existing
  `prep.rds`. Otherwise, the result of a call to
  [`prep_sdm()`](https://acanthiza.github.io/envSDM/reference/prep_sdm.md)
  with return_val = "object".

- full_run:

  Character or named list. If character, the path to an existing
  `full_run.rds`. Otherwise, the result of a call to
  [`run_full_sdm()`](https://acanthiza.github.io/envSDM/reference/run_full_sdm.md)
  with return_val = "object".

- out_dir:

  Character. Name of directory into which `.tif`s will be saved. Will be
  created if it does not exist.

- file_name:

  Character. Name to give the output prediction .tif.

- use_env_naming:

  Logical. If `TRUE`, and `is_env_pred` is `TRUE`, naming will ignore
  `file_name` and instead generate a name matching `name_env_tif()` with
  `layer` being `this_taxa` from `prep` and `start_date` being the
  minimum available `start_date` from the predictors. `pred` appears
  between `this_taxa` and `start_date`.

- predictors:

  Character. Vector of paths to predictor `.tif` files.

- is_env_pred:

  Logical. Does the naming of the directory and files in `predictors`
  follow the pattern required by `envRaster::parse_env_tif()`?

- terra_options:

  Passed to
  [`terra::terraOptions()`](https://rspatial.github.io/terra/reference/terraOptions.html).
  e.g. list(memfrac = 0.6)

- doClamp:

  Passed to
  [`terra::predict()`](https://rspatial.github.io/terra/reference/predict.html)
  (which then passes as `...` to `fun`). Possibly orphaned from older
  envSDM?

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
  and delete them if they do. Useful after a crash during predict.

- handle_errors:

  Logical. Use purrr::safely when predicting, enabling the capture of
  (m)any errors (which are then written to the log). Suggest turning off
  (i.e. `handle_errors = FALSE`) when running in a targets pipeline.

- ...:

  Passed to `...` in
  [`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html) -
  the last step in the `envSDM::predict_sdm` process. Used to provide
  additional arguments to
  [`terra::writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.html).

## Value

Character path to predicted file, usually 'pred.tif'. Output .tif and
.log, written to `out_dir`.

## Examples

``` r
  # setup -------
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  # predictors -------
  preds <- fs::dir_ls(fs::path(out_dir, "tif"))


  # data ------
  data <- fs::path(system.file(package = "envSDM"), "examples") |>
    fs::dir_ls(regexp = "prep\\.rds$"
               , recurse = TRUE
               ) |>
    tibble::enframe(name = NULL, value = "prep") |>
    dplyr::mutate(taxa = gsub("\\.rds", "", basename(dirname(prep)))
                  , tune = gsub("prep", "tune", prep)
                  , out_dir = fs::path(out_dir, taxa, "combo")
                  )

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
                                         , force_new = FALSE
                                         )
               )
#> [[1]]
#> /home/nwilloug/temp/RtmpZmDzsp/temp_libpath19cb0544a4dd6f/envSDM/examples/chg/combo/full_run.rds
#> 
#> [[2]]
#> /home/nwilloug/temp/RtmpZmDzsp/temp_libpath19cb0544a4dd6f/envSDM/examples/mjs/combo/full_run.rds
#> 
#> [[3]]
#> /home/nwilloug/temp/RtmpZmDzsp/temp_libpath19cb0544a4dd6f/envSDM/examples/wjb/combo/full_run.rds
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
                                     , force_new = FALSE
                                     )
               )

  ## .pngs -------
  if(FALSE) {

    # not working for binary 'thresh' rasters (due to an issue with number of 'classes"?):
      # Error in if (any(na.omit(x) < min(breaks)) && show.warnings) warning("Values have found that are less than the lowest break",  :
      # missing value where TRUE/FALSE needed
    purrr::walk2(data$out_dir
                 , data$out_dir
                 , \(x, y) png_from_preds(pred_dir = x
                                          , tune_dir = y
                                          , trim = FALSE
                                          , recurse = 1
                                          )
                 )

  }

  ## visualise-------
  ### mask -------
  purrr::walk(data$out_dir[file.exists(fs::path(data$out_dir, "pred.tif"))]
              , \(x) fs::path(x, "pred.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )



  # Best auc--------
  ## run full SDM --------
  data <- data %>%
    dplyr::mutate(out_dir = gsub("combo", "auc_po", out_dir))

  purrr::pmap(list(data$prep
                    , data$tune
                    , data$out_dir
                    )
               , \(a, b, c) run_full_sdm(prep = a
                                         , tune = b
                                         , out_dir = c
                                         , use_metric = "auc_po"

                                         # passed to tune_sdm via dots
                                         , metrics_df = envSDM::sdm_metrics
                                         , force_new = FALSE
                                         )
              )
#> [[1]]
#> /home/nwilloug/temp/RtmpZmDzsp/temp_libpath19cb0544a4dd6f/envSDM/examples/chg/auc_po/full_run.rds
#> 
#> [[2]]
#> /home/nwilloug/temp/RtmpZmDzsp/temp_libpath19cb0544a4dd6f/envSDM/examples/mjs/auc_po/full_run.rds
#> 
#> [[3]]
#> /home/nwilloug/temp/RtmpZmDzsp/temp_libpath19cb0544a4dd6f/envSDM/examples/wjb/auc_po/full_run.rds
#> 

  ## predict -------
  purrr::pwalk(list(data$prep
                    , data$out_dir
                    )
               , \(a, b) predict_sdm(prep = a
                                     , full_run = fs::path(b, "full_run.rds")
                                     , out_dir = b
                                     , predictors = preds
                                     , is_env_pred = FALSE
                                     , check_tifs = TRUE
                                     , force_new = FALSE
                                     )
               )


  ## visualise-------
  ### mask -------
  purrr::walk(data$out_dir[file.exists(fs::path(data$out_dir, "pred.tif"))]
              , \(x) fs::path(x, "pred.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )





```
