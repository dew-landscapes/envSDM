
#' Threshold a previously predicted SDM
#'
#'
#' @param pred_file Character. File path of predicted sdm to threshold.
#' @param threshold Numeric. > 0 and < 1. Threshold to apply to the raster
#' stored in the file at `pred_file`. Often this value will be available within
#' the result of a call to `tune_sdm()`. e.g. `mod <- rio::import("tune.rds")`
#' and then `mod$e[[1]]@thresholds$max_spec_sens`
#' @param out_dir Character. Name of directory into which `.tif`s will be saved.
#' Will be created if it does not exist. Defaults to `dirname(pred_file)`
#' @param thresh_file Character. Name to give the output threshold. Defaults to
#' `gsub("pred", "thresh", pred_file)`
#' @param terra_options Passed to `terra::terraOptions()`. e.g. list(memfrac =
#' 0.6)
#' @param force_new Logical. If output files already exist, should they be
#' remade?
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful to keep RAM use down when running SDMs for many, many taxa,
#' especially if done in parallel.
#' @param check_tifs Logical. Check if any output `.tif` files error on
#' `terra::rast()` and delete them if they do. Useful after a crash during
#' pred_file.
#'
#' @return List. `list(thresh = thresh_file)` and corresponding file written.
#'
#' @export
#'
#' @example inst/examples/thresh_sdm_ex.R
#'
thresh_sdm <- function(pred_file
                       , threshold
                       , out_dir = dirname(pred_file)
                       , thresh_file = gsub("__pred__", "__thresh__", pred_file)
                       , terra_options
                       , force_new
                       , do_gc
                       , check_tifs
                       ) {

  # setup -------
  ## start timer ------
  start_thresh <- Sys.time()

  ## this taxa ------
  this_taxa <- stringr::str_extract(basename(pred_file), "[a-zA-Z]+\\s[a-zA-Z]+")

  ## files -----
  ### new -------
  thresh_file <- fs::path(out_dir, thresh_file)

  log_file <- gsub("tif$", "log", pred_file) |>
      gsub("__pred__", "__thresh_log__", x = _)

  ### out_dir ------
  fs::dir_create(dirname(thresh_file))

  if(!dir.exists(dirname(thresh_file))) {

    stop("can't create ", dirname(thresh_file))

  }

  ### log --------
  readr::write_lines(paste0("\n"
                            , this_taxa
                            , "\threshold started at "
                            , start_time
                            )
                     , file = log_file
                     , append = TRUE
                     )

  ## test tifs ------
  # need to test before running, in case the test deletes an incomplete .tif
  if(check_tifs) {

    safe_rast <- purrr::safely(terra::rast)

    tests <- if(file.exists(thresh_file)) safe_rast(thresh_file) else NULL

    if(length(tests)) {

      warning(thresh_file
              , " will be deleted as it errored on terra::rast()"
              )

      fs::file_delete(thresh_file)

    }

  }

  # run?-----
  run <- if(file.exists(thresh_file)) force_new else TRUE # output not already exists (unless force_new)

  if(run) {

    safe_app <- purrr::safely(terra::)

    t <- safe_app(terra::rast(pred_file)
                  , \(i) i > thresh
                  , filename = thresh_file
                  , overwrite = TRUE
                  , wopt = list(datatype = "INT1U"
                                , names = this_taxa
                                )
                  )

    if(is.null(t$error)) {

      readr::write_lines(paste0("thresh done in "
                                , round(difftime(Sys.time(), start_thresh, units = "mins"), 2)
                                , " minutes"
                                )
                         , file = log_file
                         , append = TRUE
                         )

      rm(t)

      if(do_gc) {

        gc()

      }

      res <- list(thresh = thresh_file)

    } else {

      m <- as.character(t$error)

      readr::write_lines(paste0(m
                                , "\n"
                                , "threshold file not completed"
                                )
                         , file = log_file
                         , append = TRUE
                         )

      stop(m)

    }

  } else {

    message("threshold file: "
            , thresh_file
            , " already exists"
            )

    res <- list(thresh = thresh_file)

  }

  return(res)

}







