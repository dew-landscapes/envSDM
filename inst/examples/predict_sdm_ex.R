
  # setup -------
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  env_dat <- system.file("ex/bio.tif", package = "predicts")


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
                                         )
               )


  ## predict -------
  purrr::pwalk(list(data$prep
                    , data$out_dir
                    )
               , \(a, b) predict_sdm(prep = a
                                     , full_run = fs::path(b, "full_run.rds")
                                     , out_dir = b
                                     , predictors = env_dat
                                     , is_env_pred = FALSE
                                     , limit_to_boundary = TRUE
                                     , check_tifs = TRUE
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
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "mask.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )

  ### threshold ------
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "thresh.tif") %>%
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
                                         )
                                     )
                 )


  ## predict -------
  purrr::pwalk(list(data$prep
                    , data$out_dir
                    )
               , \(a, b) predict_sdm(prep = a
                                     , full_run = fs::path(b, "full_run.rds")
                                     , out_dir = b
                                     , predictors = env_dat
                                     , is_env_pred = FALSE
                                     , limit_to_boundary = TRUE
                                     , check_tifs = TRUE
                                     )
               )


  ## visualise-------
  ### mask -------
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "mask.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )

  ### threshold ------
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "thresh.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )




