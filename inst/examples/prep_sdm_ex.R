
out_dir <- file.path(system.file(package = "envSDM"), "examples")

# data ---------
hold_prop <- c(0, 0.3)
stretch <- c(10, 30)
new_bg_test <- c(T, F)

sdms <- expand.grid(hold_prop = hold_prop, stretch = stretch, new_bg_test = new_bg_test)

data <- fs::dir_ls(out_dir, regexp = "\\.csv$")[[1]] |>
  tibble::enframe(name = NULL, value = "path") |>
  dplyr::mutate(taxa = gsub("\\.csv", "", basename(path))) |>
  dplyr::cross_join(sdms) |>
  tidyr::unite(col = "out_dir"
               , tidyselect::any_of(names(sdms))
               , sep = "__"
               , remove = FALSE
               ) |>
  dplyr::mutate(out_dir = fs::path(dirname(path), out_dir)
                , out_mcp = fs::path(out_dir, "mcp.parquet")
                , dens_ras = fs::path(out_dir, "density.tif")
                , prep = fs::path(out_dir, "prep.rds")
                , tune = fs::path(out_dir, "tune.rds")
                , full_run = fs::path(out_dir, "full_run.rds")
                , pred = fs::path(out_dir, "pred.tif")
                , thresh = fs::path(out_dir, "thresh.tif")
                )

# save here for use in other examples
saveRDS(data, fs::path(out_dir, "data.rds"))

# predictors -------
preds <- fs::dir_ls(fs::path(out_dir, "tif"))

# clip --------
# make a clip boundary so mcps stay terrestrial
clip <- terra::as.polygons(terra::rast(preds)[[1]] > -Inf) |>
  sf::st_as_sf()

# mcps --------

purrr::pwalk(list(data$path
                 , data$out_mcp
                 )
             , \(x, y) envDistribution::make_mcp(readr::read_csv(x), y, pres_x = "cell_long", pres_y = "cell_lat"
                                                 , clip = clip
                                                 , dens_int = 50000
                                                 )
             )


# prep -----------
# use the just created mcps (this allows using, say, a different spatial reliability threshold for the mcps)

purrr::pwalk(list(data$taxa
                  , data$out_dir
                  , data$path
                  , data$out_mcp
                  , data$hold_prop
                  , data$stretch
                  )
             , function(a, b, c, d, e, f) prep_sdm(this_taxa = a
                                             , out_dir = b
                                             , presence = readr::read_csv(c)
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
                                             , num_bg = 100 # way too few
                                             , bg_prop_cells = 0.01 # but ensure there are at least 1% of cells with backgrounds
                                             #, force_new = TRUE
                                             )
             )

# example of 'prep'
prep <- rio::import(data$prep[[1]], trust = TRUE)

names(prep)

# variables to remove prior to SDM
prep$reduce_env$remove


# Background points
if(require("tmap")) {

  preps <- data$prep |>
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
            , size = 0.05
            ) +
    tm_facets(by = "stretch")


}
