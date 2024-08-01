
#' Generate .png (map) files for each prediction
#'
#' Finds any .tif files in `dir` and writes them to .png files. Includes the
#' retrieval and addition to the map of: various SDM metrics; and the original
#' presence points.
#'
#' @param dir Character. Name of directory containing predicted .tif(s) to
#' save as .png
#' @param trim Logical. Trim NA values from outside (using `terra::trim()`)
#' @param force_new Logical. If .png file already exists, recreate it?
#' @param include_blocks Logical. Make a .png (map) of the blocks?
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param ... Passed to `fs::dir_ls()`
#'
#' @return `invisible(NULL)`. Writes .png files with the same file name as any
#' .tif files
#' @export
#'
#' @examples
  png_from_preds <- function(dir
                              , trim = TRUE
                              , force_new = FALSE
                              , include_blocks = TRUE
                              , do_gc = TRUE
                              , ...
                              ) {

    tifs <- fs::dir_ls(dir
                       , regexp = "tif$"
                       , ...
                       )

    if(length(tifs)) {

      tif_dirs <- tibble::enframe(unique(dirname(tifs))
                                  , name = NULL
                                  , value = "path"
                                  ) %>%
        dplyr::mutate(meta = gsub(paste0(dir,"|\\/"), "", path))

      pred_tifs <- tif_dirs %>%
        dplyr::filter(meta != "") %>%
        dplyr::mutate(tune = purrr::map(path, \(x) rio::import(fs::path(x, "tune.rds")))) %>%
        tidyr::unnest(cols = c(tune)) %>%
        dplyr::mutate(auc = purrr::map_dbl(e, \(x) x@stats$auc_po)
                      , CBI_r = purrr::map_dbl(e, \(x) x@stats$CBI_rescale)
                      , IMAE = purrr::map_dbl(e, \(x) x@stats$IMAE)
                      , combo = auc * CBI_r * IMAE
                      , tr = purrr::map_dbl(e, \(x) x@thresholds$max_spec_sens)
                      ) %>%
        dplyr::mutate(tif = purrr::map(path, \(x) fs::dir_ls(x, regexp = "tif$"))) %>%
        tidyr::unnest(cols = c(tif)) %>%
        dplyr::mutate(png = gsub("tif$", "png", tif)
                      , todo = any(!file.exists(png), force_new)
                      )

      if(any(pred_tifs$todo)) {

        pred_png <- function(tif
                             , tune_args
                             , auc
                             , CBI_r
                             , IMAE
                             , combo
                             , tr
                             , png
                             , recs
                             ) {

          title <- paste0(basename(dirname(tif))
                          , "."
                          , tune_args
                          , paste0(".a:"
                                   , round(auc, 2)
                                   )
                          , paste0(".c:"
                                   , round(CBI_r, 2)
                                   )
                          , paste0(".i:"
                                   , round(IMAE, 2)
                                   )
                          , paste0(".tr:"
                                   , round(tr, 2)
                                   )
                          )

          png(filename = png
              , width = 8
              , height = 8
              , units = "in"
              , bg = "white"
              , res = 300
              )

          terra::plot(if(trim) terra::trim(terra::rast(tif)) else terra::rast(tif)
                      , main = title
                      , fun = function() {

                        terra::points(recs
                                      , pch = 4
                                      , cex = 0.5
                                      , col = scales::alpha("blue", 0.5)
                                      )

                        }
                      )

          dev.off()

        }

        recs <- rio::import(fs::path(dir, "prep.rds"))$presence %>%
          sf::st_as_sf(coords = c("x", "y")
                       , crs = terra::crs(terra::rast(tifs[[1]]))
                       ) %>%
          sf::st_coordinates()

        purrr::pwalk(list(pred_tifs$tif[pred_tifs$todo]
                          , pred_tifs$tune_args[pred_tifs$todo]
                          , pred_tifs$auc[pred_tifs$todo]
                          , pred_tifs$CBI_r[pred_tifs$todo]
                          , pred_tifs$IMAE[pred_tifs$todo]
                          , pred_tifs$combo[pred_tifs$todo]
                          , pred_tifs$tr[pred_tifs$todo]
                          , pred_tifs$png[pred_tifs$todo]
                          )
                     , pred_png
                     , recs = recs
                     )




      }

    }

    if(all(length(tifs), include_blocks)) {

      block_file <- fs::path(dir, "blocks.png")

      if(any(!file.exists(block_file), force_new)) {

        prep_file <- fs::path(dir, "prep.rds")

        if(any(exists("prep"), file.exists(prep_file))) {

          if(!exists("prep")) {

            prep <- rio::import(prep_file)

          }

          title <- paste0(basename(dir)
                          , "\nSpatial folds = "
                          , prep$spatial_folds_used
                          , "\nFolds = "
                          , length(unique(prep$blocks$block))
                          , "\nPresences = "
                          , format(nrow(prep$blocks[prep$blocks$pa == 1,]), big.mark = ",")
                          , "\nBackground = "
                          , format(nrow(prep$blocks[prep$blocks$pa == 0,]), big.mark = ",")
                          )

          m <- prep$blocks %>%
            dplyr::select(x, y, block) %>%
            sf::st_as_sf(coords = c("x", "y")
                         , crs = sf::st_crs(terra::rast(tifs[[1]]))
                         ) %>%
            tmap::tm_shape() +
              tmap::tm_dots(col = "block"
                            , size = 0.01
                            , palette = "viridis"
                            , breaks = seq(1, max(prep$blocks$block) + 1, 1)
                            , labels = as.character(seq(1, max(prep$blocks$block), 1))
                            , alpha = 0.5
                            ) +
            tmap::tm_credits(text = title
                             , position = c("left", "bottom")
                             )

          tmap::tmap_save(m
                          , block_file
                          )

        }

      }

    }

    if(do_gc) {

      rm(list = ls())

      gc()

    }

    return(invisible(NULL))

  }
