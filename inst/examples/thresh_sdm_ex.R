
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

  ## visualise-------
  ### threshold -------
  purrr::walk(data$out_dir
              , \(x) fs::path(x, "thresh.tif") %>%
                terra::rast() %>%
                terra::trim() %>%
                terra::plot()
              )






