
  library("envSDM")

  out_dir <- file.path(system.file(package = "envSDM"), "examples")

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

  env_dat <- system.file("ex/bio.tif", package = "predicts")

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
                                            , dens_res = 1000 # ignored as decimal degrees preds
                                            , use_ecdf = TRUE
                                            )
               )

  prep_ex <- rio::import(fs::path(data$out_dir[[1]], "prep.rds"))

  names(prep_ex)

  # Predict boundary with presences on 'land'
  m <- terra::plot(!is.na(terra::rast(env_dat)[[1]]))

  m <- m %>%
    terra::plot(terra::vect(prep_ex$predict_boundary)
                , add = TRUE
                , col = 5
                , alpha = 0.5
                )

  presences <- prep_ex$presence %>%
    sf::st_as_sf(coords = c("x", "y")
                 , crs = sf::st_crs(terra::rast(env_dat[[1]]))
                 ) %>%
    terra::vect()

  m <- m %>%
    terra::plot(presences, add = TRUE)

  m

  # Spatial blocks
  head(prep_ex$blocks)

  background <- prep_ex$blocks %>%
    dplyr::filter(pa == 0) %>%
    sf::st_as_sf(coords = c("x", "y")
                 , crs = sf::st_crs(terra::rast(env_dat[[1]]))
                 ) %>%
    terra::vect()

  terra::plot(background
              , pch = 4
              , cex = 0.5
              , col = scales::alpha("blue", 0.5)
              , add = TRUE
              )
