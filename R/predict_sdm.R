
#' Predict from SDM
#'
#' The resulting `pred.tif` is masked to the boundary provided to the
#' `pred_limit` argument of prep_sdm; or generated in prep_sdm from the
#' `pred_limit`, `limit_buffer` and `pred_clip` arguments.
#'
#' @param prep Character or named list. If character, the path to an existing
#' `prep.rds`. Otherwise, the result of a call to `prep_sdm()` with return_val =
#' "object".
#' @param full_run Character or named list. If character, the path to an
#' existing `full_run.rds`. Otherwise, the result of a call to `run_full_sdm()`
#' with return_val = "object".
#' @param out_dir Character. Name of directory into which `.tif`s will be saved.
#' Will be created if it does not exist.
#' @param file_name Character. Name to give the output prediction .tif.
#' @param use_env_naming Logical. If `TRUE`, and `is_env_pred` is `TRUE`, naming
#' will ignore `file_name` and instead generate a name matching
#' `name_env_tif()` with `layer` being `this_taxa` from `prep` and `start_date`
#' being the minimum available `start_date` from the predictors. `pred` appears
#' between `this_taxa` and `start_date`.
#' @param predictors Character. Vector of paths to predictor `.tif` files.
#' @param is_env_pred Logical. Does the naming of the directory and files in
#' `predictors` follow the pattern required by `envRaster::parse_env_tif()`?
#' @param terra_options Passed to `terra::terraOptions()`. e.g. list(memfrac =
#' 0.6)
#' @param doClamp Passed to `terra::predict()` (which then passes as `...` to
#' `fun`). Possibly orphaned from older envSDM?
#' @param force_new Logical. If output files already exist, should they be
#' remade?
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful to keep RAM use down when running SDMs for many, many taxa,
#' especially if done in parallel.
#' @param check_tifs Logical. Check if any output `.tif` files error on
#' `terra::rast()` and delete them if they do. Useful after a crash during
#' predict.
#' @param handle_errors Logical. Use purrr::safely when predicting, enabling the
#' capture of (m)any errors (which are then written to the log). Suggest turning
#' off (i.e. `handle_errors = FALSE`) when running in a targets pipeline.
#' @param ... Passed to `...` in `terra::mask()` - the last step in the
#' `envSDM::predict_sdm` process. Used to provide additional arguments to
#' `terra::writeRaster`.
#'
#' @return Character path to predicted file, usually 'pred.tif'. Output .tif
#' and .log, written to `out_dir`.
#'
#' @export
#'
#' @example inst/examples/predict_sdm_ex.R
#'
predict_sdm <- function(prep
                        , full_run
                        , out_dir
                        , file_name = "pred.tif"
                        , use_env_naming = FALSE
                        , predictors = NULL
                        , is_env_pred = FALSE
                        , terra_options = NULL
                        , doClamp = TRUE
                        , force_new = FALSE
                        , do_gc = FALSE
                        , check_tifs = FALSE
                        , handle_errors = TRUE
                        , ...
                        ) {

  # setup -------
  ## start timer ------
  start_time <- Sys.time()

  ## prep -------
  if(! "list" %in% class(prep)) prep <- rio::import(prep, trust = TRUE)

  ## tune ---------
  if(! "list" %in% class(full_run)) full_run <- rio::import(full_run, trust = TRUE)

  ## this taxa ------
  this_taxa <- prep$this_taxa

  ## files -----

  ### predictors -------
  if(is_env_pred) {

    pred_df <- envRaster::name_env_tif(tibble::tibble(path = predictors), parse = TRUE)

    pred_names <- pred_df %>%
      dplyr::pull(name)

    min_date <- min(as.Date(pred_df$start_date))

    x <- terra::rast(predictors)

    names(x) <- pred_names

    if(use_env_naming) {

      file_name <- fs::path(unique(pred_df$period)
                            , paste0(this_taxa, "__", gsub("\\.tif", "", file_name), "__", min_date, ".tif")
                            )

    }

  } else {

    x <- terra::rast(predictors)

  }

  ### new -------
  pred_file <- fs::path(out_dir, file_name[1])
  log_file <- gsub("\\.tif$", ".log", pred_file)

  ### out_dir ------
  fs::dir_create(dirname(pred_file))

  if(!dir.exists(dirname(pred_file))) {

    stop("can't create ", dirname(pred_file))

  }

  ## test tifs ------
  # need to test before running, in case the test deletes an incomplete .tif
  if(check_tifs) {

    safe_rast <- purrr::safely(terra::rast)

    tests <- if(file.exists(pred_file)) safe_rast(pred_file)

    if(!is.null(tests$error)) {

      warning(pred_file
              , " will be deleted as it errored on terra::rast()"
              )

      fs::file_delete(pred_file)

    }

  }

  # run?-----
  run <- all(!prep$abandoned
             , prep$finished
             , full_run$finished
             , if(file.exists(pred_file)) force_new else TRUE # output not already exists (unless force_new)
             )

  if(run) {

    ### log --------
    log <- paste0(this_taxa
                  , "\npredict process started at "
                  , start_time
                  )

    algo <- full_run$tune_mean$algo[[1]]

    mod <- full_run[paste0("tune_", algo)][[1]]

    if(nrow(mod) > 1) {

      note <- paste0(nrow(mod)
                     , " tunes present. Using the tune from the first row."
                     )

      warning(note)

      log <- paste0(log
                    , "\n"
                    , note
                    )

    }

    # pred -------

    run <- if(file.exists(pred_file)) force_new else TRUE

    if(all(!is.null(predictors), run)) {

      gc()

      if(algo == "maxnet") requireNamespace("maxnet", quietly = TRUE)
      if(algo == "rf") requireNamespace("randomForest", quietly = TRUE)
      if(algo == "env") requireNamespace("predicts", quietly = TRUE)

      ## terra options -------
      if(!is.null(terra_options)) {

        do.call(terra::terraOptions
                , args = terra_options
                )

      }

      ## predict_stack---------
      use_boundary <- if(! identical(terra::crs(x), terra::crs(prep$predict_boundary))) {

        prep$predict_boundary |>
          terra::vect() |>
          terra::densify(50000) |>
          sf::st_as_sf() |>
          sf::st_transform(crs = sf::st_crs(x)) |>
          sf::st_make_valid()

      } else prep$predict_boundary

      terra::window(x) <- terra::ext(terra::vect(use_boundary))

      ## predict--------
      pred_start <-  Sys.time()

      safe_predict <- purrr::safely(terra::predict)

      m <- paste0("predict for "
                  , this_taxa
                  , " from "
                  , algo
                  , " model with arguments: "
                  , gsub("\\.", ",", gsub(":", " =", full_run$tune_mean$tune_args))
                  , "\n out file will be "
                  , basename(pred_file)
                  )

      message(m)

      log <- paste0(log
                    , "\n"
                    , m
                    )

      window_predict_file <- paste0(tempfile(), ".tif")

      p <- safe_predict(object = x
                        , model = mod$m[[1]]
                        , type = if(algo == "rf") "prob" else "cloglog"
                        , clamp = doClamp
                        , na.rm = TRUE
                        , index = if(algo == "rf") "1" else NULL
                        , names = this_taxa
                        , filename = window_predict_file
                        , overwrite = TRUE
                        )

      ## if error -------
      if(!is.null(p$error)) {

        m <- paste0("error: "
                    , as.character(p$error)
                    , ". Predict process not finished."
                    )

        message(m)

        log <- paste0(log
                      , "\n"
                      , m
                      )

        ### handle errors -------
        if(! handle_errors) stop(m)

      }

      ## else mask -------
      if(is.null(p$error)) {

        log <- paste0(log
                      , "\n"
                      , "predict finished in "
                      , round(difftime(Sys.time(), pred_start, units = "mins"), 2)
                      , " minutes"
                      )

        predict_mask_start <- Sys.time()

        m <- paste0("predict mask for ", this_taxa)

        message(m)

        log <- paste0(log
                      , "\n"
                      , m
                      )

        terra::mask(p$result
                    , mask = terra::vect(use_boundary)
                    , filename = pred_file
                    , overwrite = TRUE
                    , ...
                    )

        if(file.exists(window_predict_file)) fs::file_delete(window_predict_file)

        log <- paste0(log
                      , "\n"
                      , "predict mask created in "
                      , round(difftime(Sys.time(), predict_mask_start, units = "mins"), 2)
                      , " minutes"
                      )

        ## finish ----
        log <- paste0(log
                      , "\n"
                      , "predict process finished. elapsed time: "
                      , round(difftime(Sys.time(), start_time, units = "mins"), 2)
                      , " minutes"
                      )

      }

      if(do_gc) {

        rm(p, x)

        gc()

      }

    }

  } else {

    # if prep abandoned but pred_file exists rename it, then delete it
    # to the same for any thresh_ file
    # this can happen occasionally when adjusting the incoming data so that a taxa that previously had an SDM is now abandoned
    if(all(prep$abandoned, file.exists(pred_file))) {

      fs::file_copy(pred_file
                    , gsub(basename(pred_file), paste0("renamed__", basename(pred_file)), x = pred_file)
                    )

      fs::file_delete(pred_file)

    }

    # if pred_file does not exist, then create a log

    if(! file.exists(pred_file)) {

      log <- paste0(this_taxa
                    , "\npredict process started at "
                    , start_time
                    )

      m <- paste0("predict did not run as:\n "
                  , if(prep$abandoned) "prep abandoned\n "
                  , if(! prep$finished) "prep did not finish\n "
                  , if(any(prep$abandoned, ! prep$finished)) paste0("prep log:\n  ", paste0(prep$log, collapse = "\n  "), "\n ")
                  , if (! full_run$finished) paste0("full run did not finish\n full run log:\n  ", paste0(full_run$log, collapse = "\n  "), "\n ")
                  , "\n"
                  )

      log <- paste0(log
                    , "\n"
                    , m
                    )

    }

  }

  # write log ------
  # only if it exists
  if(exists("log", inherits = FALSE)) {

    writeLines(log
               , log_file
               )

  }

  if(do_gc) {

    stuff <- ls()

    delete_stuff <- stuff[!grepl("_file$", stuff)]

    rm(list = delete_stuff)

    gc()

  }

  res <- if(file.exists(pred_file)) pred_file else NULL

  return(res)

}

