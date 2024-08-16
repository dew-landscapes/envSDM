
#' Generate .png (map) files for each prediction
#'
#' Finds any .tif files in `pred_dir` and writes them to .png files. Includes the
#' retrieval and addition to the map of: various SDM metrics; and the original
#' presence points.
#'
#' @param pred_dir Character. Name of directory containing predicted .tif(s) to
#' save as .png
#' @param tune_dir Character. Name of directory containing the tune results used
#' to make the predictions in pred_dir.
#' @param trim Logical. Trim NA values from outside (using `terra::trim()`)
#' @param force_new Logical. If .png file already exists, recreate it?
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param ... Passed to `fs::dir_ls()`
#'
#' @return `invisible(NULL)`. Writes .png files with the same file name as any
#' .tif files
#' @export
#' @example inst/examples/predict_sdm_ex.R
png_from_preds <- function(pred_dir
                           , tune_dir = NULL
                           , trim = TRUE
                           , force_new = FALSE
                           , do_gc = TRUE
                           , ...
                           ) {

  if(is.null(tune_dir)) tune_dir = pred_dir

  tifs <- fs::dir_ls(pred_dir
                     , regexp = "tif$"
                     , ...
                     )

  if(length(tifs)) {

    tif_dirs <- tibble::enframe(unique(dirname(tifs))
                                , name = NULL
                                , value = "pred_dir"
                                ) %>%
      dplyr::mutate(metric = basename(pred_dir)
                    , tune_dir = fs::path(tune_dir, metric)
                    ) %>%
      dplyr::filter(metric %in% c("combo", envSDM::sdm_metrics$metric)) %>%
      dplyr::mutate(taxa = basename(dirname(pred_dir)))

    pred_tifs <- tif_dirs %>%
      dplyr::filter(metric != "") %>%
      dplyr::mutate(tune_dir = fs::path(tune_dir)
                    , tune = purrr::map(tune_dir
                                      , \(x) rio::import(fs::path(x, "tune.rds"))
                                      )
                    ) %>%
      tidyr::unnest(cols = c(tune)) %>%
      dplyr::mutate(auc = purrr::map_dbl(e, \(x) x@stats$auc_po)
                    , CBI_r = purrr::map_dbl(e, \(x) x@stats$CBI_rescale)
                    , IMAE = purrr::map_dbl(e, \(x) x@stats$IMAE)
                    , combo = auc * CBI_r * IMAE
                    , tr = purrr::map_dbl(e, \(x) x@thresholds$max_spec_sens)
                    ) %>%
      dplyr::mutate(tif = purrr::map(pred_dir, \(x) fs::dir_ls(x, regexp = "tif$"))) %>%
      dplyr::select(! dplyr::where(is.list), tif) %>%
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

        text <- paste0("'best' "
                       , basename(dirname(tif))
                        , "\n"
                        , paste0("tune args: "
                                 , tune_args
                                 , ". "
                                 )
                        , paste0("auc_po: "
                                 , round(auc, 2)
                                 , "\n"
                                 )
                        , paste0("CBI_rescale: "
                                 , round(CBI_r, 2)
                                 , ". "
                                 )
                        , paste0("IMAE :"
                                 , round(IMAE, 2)
                                 , ". "
                                 )
                        , paste0("Threshold: "
                                 , round(tr, 2)
                                 , "."
                                 )
                        )

        ras <- if(trim) terra::trim(terra::rast(tif)) else terra::rast(tif)
        title <- basename(dirname(dirname(terra::sources(ras))))

        m <- tmap::tm_shape(ras) +
          tmap::tm_raster(title = title) +
          tmap::tm_shape(recs) +
          tmap::tm_dots(alpha = 0.5) +
          tmap::tm_credits(text = text
                           , position = c("left", "bottom")
                           ) +
          tmap::tm_compass() +
          tmap::tm_scale_bar()

        tmap::tmap_save(m
                        , png
                        )

      }

      recs <- rio::import(fs::path(tune_dir, "prep.rds"))$presence %>%
        sf::st_as_sf(coords = c("x", "y")
                     , crs = terra::crs(terra::rast(tifs[[1]]))
                     )

      purrr::pwalk(list(tif = pred_tifs$tif[pred_tifs$todo]
                        , tune_args = pred_tifs$tune_args[pred_tifs$todo]
                        , auc = pred_tifs$auc[pred_tifs$todo]
                        , CBI_r = pred_tifs$CBI_r[pred_tifs$todo]
                        , IMAE = pred_tifs$IMAE[pred_tifs$todo]
                        , combo = pred_tifs$combo[pred_tifs$todo]
                        , tr = pred_tifs$tr[pred_tifs$todo]
                        , png = pred_tifs$png[pred_tifs$todo]
                        )
                   , pred_png
                   , recs = recs
                   )

    }

  }

  if(do_gc) {

    rm(list = ls())

    gc()

  }

  return(invisible(NULL))

}
