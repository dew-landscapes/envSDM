
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  # data ---------
  data <- fs::dir_ls(out_dir, regexp = "\\.csv$")[[1]] |>
    tibble::enframe(name = NULL, value = "path") |>
    dplyr::mutate(taxa = gsub("\\.csv", "", basename(path))
                  , presence = purrr::map(path, rio::import, setclass = "tibble", trust = TRUE)
                  ) |>
    dplyr::cross_join(tibble::tibble(hold_prop = c(0))) |>
    dplyr::cross_join(tibble::tibble(repeats = c(5))) |>
    dplyr::cross_join(tibble::tibble(stretch = c(10, 20, 100))) |>
    dplyr::mutate(taxa_dir = fs::path(out_dir, paste0(taxa, "__", hold_prop, "__", repeats, "__", stretch))
                  , out_mcp = fs::path(taxa_dir, "mcp.parquet")
                  , dens_ras = fs::path(taxa_dir, "density.tif")
                  , prep_file = fs::path(taxa_dir, "prep.rds")
                  )

  # predictors -------
  preds <- fs::dir_ls(fs::path(out_dir, "tif"))

  # clip --------
  # make a clip boundary so mcps stay terrestrial
  clip <- terra::as.polygons(terra::rast(preds)[[1]] > -Inf) |>
    sf::st_as_sf()

  # mcps --------

  purrr::pwalk(list(data$presence
                   , data$out_mcp
                   )
               , \(x, y) envDistribution::make_mcp(x, y, pres_x = "cell_long", pres_y = "cell_lat"
                                                   , clip = clip
                                                   , dens_int = 50000
                                                   )
               )


  # prep -----------
  # use the just created mcps (this allows using, say, a different spatial reliability threshold for the mcps)

  purrr::pwalk(list(data$taxa
                    , data$taxa_dir
                    , data$presence
                    , data$out_mcp
                    , data$hold_prop
                    , data$stretch
                    )
               , function(a, b, c, d, e, f) prep_sdm(this_taxa = a
                                               , out_dir = b
                                               , presence = c
                                               , pres_x = "cell_long"
                                               , pres_y = "cell_lat"
                                               , predictors = preds
                                               , is_env_pred = FALSE
                                               , pred_limit = d
                                               , limit_buffer = 10000
                                               , folds = 5
                                               , repeats = 5
                                               , hold_prop = e
                                               , dens_res = 1000 # ignored as decimal degrees preds
                                               , reduce_env_thresh_corr = 0.95
                                               , reduce_env_quant_rf_imp = 0.2
                                               , stretch_value = f
                                               , num_bg = 1000
                                               #, force_new = TRUE
                                               )
               )

  # example of 'prep'
  prep <- rio::import(data$prep_file[[1]], trust = TRUE)

  names(prep)

  # variables to remove prior to SDM
  prep$reduce_env$remove


  # Background points
  if(require("tmap")) {

    preps <- data$prep_file |>
      purrr::map(readRDS)

    b <- preps |>
      purrr::map("bg_points") |>
      purrr::set_names(data$stretch) |>
      dplyr::bind_rows(.id = "stretch") |>
      dplyr::mutate(stretch = as.numeric(stretch)
                    , type = "background"
                    )

    p <- preps |>
      purrr::map(\(x) sf::st_as_sf(x$presence_ras
                                   , coords = c("x", "y")
                                   , crs = x$epsg_out
                                   )
                 ) |>
      purrr::set_names(data$stretch) |>
      dplyr::bind_rows(.id = "stretch") |>
      dplyr::mutate(stretch = as.numeric(stretch)
                    , type = "presence"
                    )



    tm_shape(dplyr::bind_rows(b, p)) +
      tm_dots(fill = "type"
              , size = 0.5
              ) +
      tm_facets(by = "stretch")


  }
