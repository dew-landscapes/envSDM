
#' Generate .png (map) files for each prediction
#'
#' Finds any .tif files in `pred_dir` and writes them to .png files. Includes the
#' retrieval and addition to the map of: various SDM metrics; and the original
#' presence points.
#'
#' @param prep Character or named list. If character, the path to an existing
#' `prep.rds`. Otherwise, the result of a call to `prep_sdm()` with return_val =
#' "object".
#' @param full_run Character or named list. If character, the path to an
#' existing `full_run.rds`. Otherwise, the result of a call to `run_full_sdm`()
#' with return_val = "object".
#' @param out_dir Character. Name of directory into which `.pngs`s will be
#' saved. Will be created if it does not exist.
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
                           , full_run_dir = NULL
                           , trim = TRUE
                           , force_new = FALSE
                           , do_gc = TRUE
                           , ...
                           ) {

  # setup -------
  if(isFALSE(out_dir)) out_dir <- tempfile()

  ## out_dir ------
  if(is.character(out_dir)) {

    fs::dir_create(out_dir)

    if(!dir.exists(out_dir)) {

      stop("can't create out_dir")

    }

  }

  ## prep -------
  if(! "list" %in% class(prep)) prep <- rio::import(prep)

  ## tune ---------
  if(! "list" %in% class(full_run)) full_run <- rio::import(full_run)

  if(dir.exists(pred_dir)) {

    tifs <- fs::dir_ls(pred_dir
                       , regexp = "tif$"
                       , ...
                       )

    pngs <- gsub("tif$", "png", tifs)

  } else {

    tifs <- NULL

    mes <- paste0("The directory: "
                  , pred_dir
                  , " does not exist"
                  )

  }

  if(length(tifs)) {

    tune <- full_run$tune_mean

    todo <- any(!file.exists(pngs), force_new)

    if(todo) {

      pred_png <- function(tif
                           , tune
                           , png
                           , recs
                           ) {

        text <- paste0("'best' tune arguments with "
                       , basename(dirname(tif))
                       , " metric:"
                        , "\n"
                        , paste0(tune$tune_args
                                 , ". "
                                 )
                        , paste0("auc_po: "
                                 , round(tune$auc_po, 2)
                                 , "\n"
                                 )
                        , paste0("CBI_rescale: "
                                 , round(tune$CBI_rescale, 2)
                                 , ". "
                                 )
                        , paste0("IMAE :"
                                 , round(tune$IMAE, 2)
                                 , ". "
                                 )
                        , paste0("Threshold: "
                                 , round(tune$max_spec_sens, 2)
                                 , "."
                                 )
                        )

        ras <- if(trim) terra::trim(terra::rast(tif)) else terra::rast(tif)
        title <- basename(dirname(dirname(terra::sources(ras))))

        classes <- if(grepl("thresh", tif)) 1 else 10

        m <- tmap::tm_shape(ras) +
          tmap::tm_raster(n = classes) +
          tmap::tm_shape(recs) +
          tmap::tm_dots(alpha = 0.5) +
          tmap::tm_credits(text = text
                           , bg.color = "grey"
                           , bg.alpha = 0.5
                           ) +
          tmap::tm_graticules() +
          tmap::tm_layout(title = title
                          , legend.bg.color = "grey"
                          , legend.bg.alpha = 0.2
                          )

        tmap::tmap_save(m
                        , png
                        )

      }

      recs <- prep$presence %>%
        sf::st_as_sf(coords = c("x", "y")
                     , crs = terra::crs(terra::rast(tifs[[1]]))
                     )

      purrr::pwalk(list(tifs
                        , pngs
                        )
                   , \(x, y) pred_png(tif = x
                                      , tune = tune
                                      , png = y
                                      , recs = recs
                                      )
                   )

    }

  } else {

    if(! exists("mes")) mes <- paste0("No tifs found in ", pred_dir)

    message(mes)

  }

  if(do_gc) {

    rm(list = ls())

    gc()

  }

  return(invisible(NULL))

}
