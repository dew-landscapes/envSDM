
  # setup -------
  source(fs::path("inst", "examples", "prep_sdm_ex.R")) # need 'data' object

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

  ## visualise-------

  tifs <- data$thresh

  names <- paste0("hold_prop "
                  , data$hold_prop
                  , "; stretch "
                  , data$stretch
                  , "; new_bg "
                  , data$new_bg_test
                  )

  r <- terra::rast(tifs)
  names(r) <- names
  terra::plot(r, cex.main = 0.7, nc = 2)
