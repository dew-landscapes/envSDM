
  # setup -------
  source(fs::path("inst", "examples", "prep_sdm_ex.R")) # need 'data' object

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
