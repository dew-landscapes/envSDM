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
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0/combo/full_run.rds
#> 
#> [[2]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3/combo/full_run.rds
#> 
#> [[3]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3__1/combo/full_run.rds
#> 
#> [[4]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3__5/combo/full_run.rds
#> 
#> [[5]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0__1/combo/full_run.rds
#> 
#> [[6]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0__5/combo/full_run.rds
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
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0/auc_po/full_run.rds
#> 
#> [[2]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3/auc_po/full_run.rds
#> 
#> [[3]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3__1/auc_po/full_run.rds
#> 
#> [[4]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0.3__5/auc_po/full_run.rds
#> 
#> [[5]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0__1/auc_po/full_run.rds
#> 
#> [[6]]
#> /home/nwilloug/temp/RtmprTkVGw/temp_libpath19fa091c7edf0e/envSDM/examples/chg__0__5/auc_po/full_run.rds
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
