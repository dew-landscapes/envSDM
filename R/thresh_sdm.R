
#' Threshold a previously predicted SDM
#'
#'
#' @param pred_file Character. File path of predicted sdm to threshold.
#' @param this_taxa Character. If left as default `NULL` an attempt will be made
#' to extract a taxa name from `pred_file`
#' @param threshold Numeric. > 0 and < 1. Threshold to apply to the raster
#' stored in the file at `pred_file`. Often this value will be available within
#' the result of a call to `tune_sdm()`. e.g. `mod <- rio::import("tune.rds")`
#' and then `mod$e[[1]]@thresholds$max_spec_sens`
#' @param thresh_file Character. Name to give the output threshold. If left as
#' default `NULL`, `thresh_file` will be set to
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
#' @return Character path to threshold file, usually 'thresh.tif'. Output .tif
#' and .log, written to `out_dir`.
#'
#' @export
#'
#' @example inst/examples/thresh_sdm_ex.R
#'
thresh_sdm <- function(pred_file
                       , this_taxa = NULL
                       , threshold
                       , thresh_file = NULL
                       , terra_options = NULL
                       , force_new = FALSE
                       , do_gc = FALSE
                       , check_tifs = TRUE
                       ) {

  # pred exists -----
  if(!is.null(pred_file)) {

    # setup -------
    ## NULL args -------
    if(is.null(this_taxa)) this_taxa <- stringr::str_extract(pred_file, "[[:alpha:]]+\\s[[:alpha:]]+")
    if(is.null(thresh_file)) thresh_file <- gsub("__pred__", "__thresh__", pred_file)

    ## start timer ------
    start_thresh <- Sys.time()

    ## files -----
    ### new -------
    file_name <- basename(thresh_file)
    dir_name <- if(dirname(thresh_file) == ".") dirname(pred_file) else dirname(thresh_file)
    thresh_file <- fs::path(dir_name, file_name)

    log_file <- gsub("tif$", "log", thresh_file)

    ### log --------
    log <- paste0(this_taxa
                  , " threshold started at "
                  , start_thresh
                  )

    ## test tifs ------
    # need to test before running, in case the test deletes an incomplete .tif
    if(check_tifs) {

      safe_rast <- purrr::safely(terra::rast)

      tests <- if(file.exists(thresh_file)) safe_rast(thresh_file)

      if(!is.null(tests$error)) {

        warning(thresh_file
                , " will be deleted as it errored on terra::rast()"
                )

        fs::file_delete(thresh_file)

      }

    }

    # run?-----
    run <- if(file.exists(thresh_file)) force_new else TRUE # output not already exists (unless force_new)

    if(run) {

      m <- paste0("create "
                  , basename(thresh_file)
                  , " for "
                  , this_taxa
                  , " with threshold value: "
                  , threshold
                  )

      message(m)

      log <- paste0(log
                    , "\n"
                    , m
                    )

      safe_app <- purrr::safely(terra::app)

      t <- safe_app(terra::rast(pred_file)
                    , \(i) i > threshold
                    , filename = thresh_file
                    , overwrite = TRUE
                    , wopt = list(datatype = "INT1U"
                                  , names = this_taxa
                                  )
                    )

      if(is.null(t$error)) {

        log <- paste0(log
                      , "\n"
                      , "thresh done in "
                      , round(difftime(Sys.time(), start_thresh, units = "mins"), 2)
                      , " minutes"
                      )

        rm(t)

        if(do_gc) {

          gc()

        }

        res <- thresh_file

      } else {

        m <- as.character(t$error)

        log <- paste0(log
                      , "\n"
                      , m
                      , "\n"
                      , "threshold file not completed"
                      )

        stop(m)

      }

    } else {

      message("threshold file: "
              , thresh_file
              , " already exists"
              )

      res <- thresh_file

    }

  } else {

    message("pred_file is null")

    res <- NULL

  }

  writeLines(log
             , log_file
             )

  return(res)

}







