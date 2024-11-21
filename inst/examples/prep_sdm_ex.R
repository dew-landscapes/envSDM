
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
                  , taxa_dir = fs::path(out_dir, taxa)
                  , out_mcp = fs::path(taxa_dir, "mcp.parquet")
                  )

  env_dat <- system.file("ex/bio.tif", package = "predicts")


  # mcps --------

  # make a clip boundary so mcps stay terrestrial
  clip <- terra::as.polygons(terra::rast(env_dat)[[1]] > -Inf) |>
    sf::st_as_sf()

  purrr::pwalk(list(data$presence
                   , data$out_mcp
                   )
              , \(x, y) make_mcp(x, y, pres_x = "lon"
                                 , clip = clip
                                 )
              )


  # prep -----------
  # use the just created mcps (this allows using, say, a different spatial reliability threshold for the mcps)

  purrr::pwalk(list(data$taxa
                    , data$taxa_dir
                    , data$presence
                    , data$out_mcp
                    )
               , function(a, b, c, d) prep_sdm(this_taxa = a
                                               , out_dir = b
                                               , presence = c
                                               , pres_x = "lon"
                                               , pres_y = "lat"
                                               , predictors = env_dat
                                               , is_env_pred = FALSE
                                               , pred_limit = d
                                               , limit_buffer = 10000
                                               , dens_res = 1000 # ignored as decimal degrees preds
                                               )
               )

  # example of 'prep'
  prep <- rio::import(fs::path(data$taxa_dir[[2]], "prep.rds"))

  names(prep)

  # env variables not used
  prep$reduce_env$remove

  # Density raster
  dens_ras <- terra::rast(fs::path(data$taxa_dir[[2]], "density.tif")) %>%
    terra::mask(clip) %>%
    terra::classify(matrix(c(0, NA), ncol = 2))

  if(require("tmap")) {

    m <-
      tm_shape(clip) +
      tm_borders() +
      tm_shape(dens_ras) +
      tm_raster(title = "Background point density"
                , breaks = c(0, 2, 4, 6, 8, 10)
                , drop.levels = TRUE
                , colorNA = NULL
                ) +
      tm_legend(outside = TRUE) +
      tm_compass() +
      tm_scale_bar() +
      tm_layout(main.title = paste0("Background point density for ",  prep$inputs$this_taxa))

    m

    presences <- prep$blocks %>%
      dplyr::filter(pa == 1) %>%
      sf::st_as_sf(coords = c("x", "y")
                   , crs = 4326
                   )

    m +
      tm_shape(presences) +
        tm_dots(col = "blue")

  }

  # Background points
  if(require("tmap")) {

    blocks <- prep$blocks %>%
      dplyr::mutate(blocks = factor(block)) %>% # for map
      sf::st_as_sf(coords = c("x", "y")
                   , crs = sf::st_crs(terra::rast(env_dat[[1]]))
                   )


    tm_shape(clip) +
      tm_borders() +
      tm_shape(blocks) +
      tm_dots(title = "Background points\n coloured by block"
              , col = "block"
              ) +
      tm_legend(outside = TRUE) +
      tm_compass() +
      tm_scale_bar() +
      tm_layout(main.title = paste0("Background points for ",  prep$this_taxa))


  }
