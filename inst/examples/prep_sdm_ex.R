
  library("envSDM")
  library("tmap")

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

  prep_ex <- rio::import(fs::path(data$out_dir[[2]], "prep.rds"))

  names(prep_ex)

  # Density raster
  land <- !is.na(terra::rast(env_dat)[[1]])
  dens_ras <- terra::rast(fs::path(data$out_dir[[2]], "density.tif")) * land %>%
    terra::trim()

  m <- tm_shape(land) +
    tm_raster() +
    tm_shape(dens_ras) +
    tm_raster(title = "Presence density"
              , drop.levels = TRUE
              ) +
    tm_legend(outside = TRUE) +
    tm_compass() +
    tm_scale_bar() +
    tm_layout(main.title = paste0("Prep for ",  prep_ex$inputs$this_taxa))

  m

  # Spatial blocks
  head(prep_ex$blocks)

  blocks <- prep_ex$blocks %>%
    dplyr::mutate(blocks = factor(block)) %>% # for map
    sf::st_as_sf(coords = c("x", "y")
                 , crs = sf::st_crs(terra::rast(env_dat[[1]]))
                 )

  presences <- blocks %>%
    dplyr::filter(pa == 1)

  m +
    tm_shape(blocks) +
      tm_dots(col = "block"
              , palette = "viridis"
              , title = "Blocks"
              ) +
    tm_shape(presences) +
      tm_dots(col = "red")
