
  # setup -------
  in_dir <- file.path(system.file(package = "envSDM"), "examples")

  env_dat <- system.file("ex/bio.tif", package = "predicts")


  # data ------
  data <- file.path(system.file(package = "predicts"), "ex") |>
    fs::dir_ls(regexp = "\\.csv$") |>
    tibble::enframe(name = NULL, value = "path") |>
    dplyr::mutate(taxa = gsub("\\.csv", "", basename(path))
                  , presence = purrr::map(path, rio::import, setclass = "tibble")
                  , presence = purrr::map(presence
                                          , \(x) x |>
                                            dplyr::filter(!is.na(lat)
                                                          , !is.na(lon)
                                                          )
                                          )
                  , out_dir = fs::path(in_dir, taxa)
                  )

  # Best combo--------
  ## run full SDM --------
  purrr::pwalk(list(data$out_dir)
                 , \(a) run_full_sdm(out_dir = a
                                     , metrics_df = envSDM::sdm_metrics
                                     )
                 )


  ## predict -------
  purrr::pwalk(list(data$out_dir)
               , \(a) predict_sdm(prep_dir = a
                                  , tune_dir = fs::path(a, "combo")
                                  , predictors = env_dat
                                  , is_env_pred = FALSE
                                  , limit_to_mcp = TRUE
                                  , check_tifs = TRUE
                                  , force_new = TRUE
                                  )
               )

  ## .pngs -------
  purrr::walk2(data$out_dir
               , data$out_dir
               , \(x, y) png_from_preds(pred_dir = x
                                        , tune_dir = y
                                        , trim = FALSE
                                        #, force_new = TRUE
                                        , recurse = 1
                                        )
               )

  ## visualise-------
  ### mask -------
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "combo", "mask.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )

  ### threshold ------
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "combo", "thresh.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )

  # Best auc--------
  ## run full SDM --------
  purrr::pwalk(list(data$out_dir)
                 , \(a) run_full_sdm(out_dir = a
                                     , metrics_df = envSDM::sdm_metrics
                                     , metric = "auc_po"
                                     , save_to = fs::path(a, "auc_po")
                                     )
                 )


  ## predict -------
  purrr::pwalk(list(data$out_dir)
               , \(a) predict_sdm(prep_dir = a
                                  , tune_dir = fs::path(a, "auc_po")
                                  , predictors = env_dat
                                  , is_env_pred = FALSE
                                  , limit_to_mcp = TRUE
                                  , check_tifs = TRUE
                                  )
               )


  ## visualise-------
  ### mask -------
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "auc_po", "mask.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )

  ### threshold ------
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "auc_po", "thresh.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )




