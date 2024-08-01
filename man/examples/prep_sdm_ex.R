
  out_dir <- here::here("inst", "examples")

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
                  , out_dir = fs::path(out_dir, taxa)
                  )

  env_dat <- system.file("ex/bio.tif", package="predicts")

  purrr::pwalk(list(data$taxa
                    , data$out_dir
                    , data$presence
                    )
               , function(a, b, c) prep_sdm(this_taxa = a
                                            , out_dir = b
                                            , presence = c
                                            , pres_x = "lon"
                                            , pres_y = "lat"
                                            , predictors = env_dat
                                            , is_env_pred = FALSE
                                            , pred_limit = TRUE
                                            , limit_buffer = 10000
                                            , dens_res = 1000 # ignored as decimal degress preds
                                            )
               )

  fs::dir_ls(out_dir
             , recurse = TRUE
             )

