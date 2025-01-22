
#' Predict from SDM
#'
#' The resulting `pred.tif` is masked to the boundary provided to the
#' `pred_limit` argument of prep_sdm; or generated in prep_sdm from the
#' `pred_limit`, `limit_buffer` and `pred_clip` arguments. A threshold raster
#' can also be saved (saved as 'thresh.tif') - see `apply_thresh` argument. In
#' both cases the resulting files will have the same extent, resolution and crs
#' as the predictors. The default name of the files created (`pred.tif` and
#' `thresh.tif`) can be changed via the file_name argument.
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
#' If saving a threshold file (`apply_thresh` is `TRUE`), this must be of length
#' two with the first element used for the predictions and the second for the
#' threshold.
#' @param use_env_naming Logical. If `TRUE`, and `is_env_pred` is `TRUE`, naming
#' will ignore `file_name` and instead generate a name matching
#' `name_env_tif()` with `layer` being `this_taxa` from `prep` and `start_date`
#' being the minimum available `start_date` from the predictors. `pred` and
#' `thresh` appear between `this_taxa` and `start_date`.
#' @param predictors Character. Vector of paths to predictor `.tif` files.
#' @param is_env_pred Logical. Does the naming of the directory and files in
#' `predictors` follow the pattern required by `envRaster::parse_env_tif()`?
#' @param terra_options Passed to `terra::terraOptions()`. e.g. list(memfrac =
#' 0.6)
#' @param doClamp Passed to `terra::predict()` (which then passes as `...` to
#' `fun`). Possibly orphaned from older envSDM?
#' @param apply_thresh Logical. If `TRUE`, an output raster `thresh.tif` will be
#' created using the maximum of specificity + sensitivity. The threshold value
#' can be accessed within `tune.rds` as, say, `mod <- rio::import("tune.rds")`
#' and then `mod$e[[1]]@thresholds$max_spec_sens`
#' @param force_new Logical. If output files already exist, should they be
#' remade?
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful to keep RAM use down when running SDMs for many, many taxa,
#' especially if done in parallel.
#' @param check_tifs Logical. Check if any output `.tif` files error on
#' `terra::rast()` and delete them if they do. Useful after a crash during
#' predict.
#' @param ... Passed to `...` in `terra::mask()` - the last step in the
#' `envSDM::predict_sdm` process. Used to provide additional arguments to
#' `terra::writeRaster`.
#'
#' @return Named list of created .tif files, usually 'pred.tif' and
#' 'thresh.tif'. Output .tif(s) and .log, written to `out_dir`.
#'
#' @export
#'
#' @example inst/examples/predict_sdm_ex.R
#'
  predict_sdm <- function(prep
                          , full_run
                          , out_dir
                          , file_name = c("pred.tif", "thresh.tif")
                          , use_env_naming = FALSE
                          , predictors = NULL
                          , is_env_pred = FALSE
                          , terra_options = NULL
                          , doClamp = TRUE
                          , apply_thresh = TRUE
                          , force_new = FALSE
                          , do_gc = FALSE
                          , check_tifs = TRUE
                          , ...
                          ) {

    # setup -------
    ## start timer ------
    start_time <- Sys.time()

    ## prep -------
    if(! "list" %in% class(prep)) prep <- rio::import(prep)

    ## tune ---------
    if(! "list" %in% class(full_run)) full_run <- rio::import(full_run)

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
    if(apply_thresh) thresh_file <- fs::path(out_dir, file_name[2])
    log_file <- gsub("tif$", "log", pred_file) %>%
      gsub("__pred__", "__log__", .)

    ### out_dir ------
    fs::dir_create(dirname(pred_file))

    if(!dir.exists(dirname(pred_file))) {

      stop("can't create ", dirname(pred_file))

    }

    ### log --------
    readr::write_lines(paste0("\n\n"
                              , this_taxa
                              , "\npredict start at "
                              , start_time
                              )
                       , file = log_file
                       , append = TRUE
                       )

    ## test tifs ------
    # need to test before running, in case the test deletes an incomplete .tif
    if(all(check_tifs, dir.exists(out_dir))) {

      safe_rast <- purrr::safely(terra::rast)

      test_files <- fs::dir_ls(out_dir
                               , regexp = "tif$|nc$"
                               , recurse = 1
                               )

      tests <- purrr::map(test_files
                          , safe_rast
                          ) %>%
        purrr::map("error") %>%
        purrr::compact()

      if(length(tests)) {

        warning("raster files to be deleted as they errored on terra::rast(): "
                , envFunc::vec_to_sentence(names(tests)[file.exists(names(tests))])
                )

        purrr::map(names(tests)[file.exists(names(tests))]
                   , unlink
                   , force = TRUE
                   )

      }

    }


    # run?-----
    run <- all(!prep$abandoned
               , prep$finished
               , full_run$finished
               , if(file.exists(thresh_file)) force_new else TRUE # output not already exists (unless force_new)
               )

    if(run) {

      algo <- full_run$tune_mean$algo[[1]]

      mod <- full_run[paste0("tune_", algo)][[1]]

      if(nrow(mod) > 1) {

        note <- paste0(nrow(mod)
                       , " tunes present. Using the tune from the first row."
                       )

        warning(note)

        readr::write_lines(note
                           , file = log_file
                           , append = TRUE
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
        predict_stack_start <- Sys.time()

        m <- paste0("predict window for ", this_taxa)

        message(m)

        readr::write_lines(m
                           , file = log_file
                           , append = TRUE
                           )

        window(x) <- terra::ext(terra::vect(prep$predict_boundary))

        readr::write_lines(paste0("predict window created in "
                                  , round(difftime(Sys.time(), predict_stack_start, units = "mins"), 2)
                                  , " minutes"
                                  )
                           , file = log_file
                           , append = TRUE
                           )

        ## predict--------
        pred_start <-  Sys.time()

        safe_predict <- purrr::safely(terra::predict)

        m <- paste0("predict for "
                    , this_taxa
                    , " from "
                    , algo
                    , " model with arguments: "
                    , gsub("\\.", ",", gsub(":", " =", full_run$tune_mean$tune_args))
                    , " and threshold: "
                    , round(full_run$tune_mean$max_spec_sens, 2)
                    , "\n out file will be "
                    , basename(pred_file)
                    )

        message(m)

        readr::write_lines(m
                           , file = log_file
                           , append = TRUE
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

        if(!is.null(p$error)) {

          m <- paste0("error: "
                      , as.character(p$error)
                      )

          message(m)

          readr::write_lines(m
                             , file = log_file
                             , append = TRUE
                             )

        } else {

          terra::mask(p$result
                      , mask = terra::vect(prep$predict_boundary)
                      , filename = pred_file
                      , ...
                      )

          if(file.exists(window_predict_file)) fs::file_delete(window_predict_file)

        }

        if(do_gc) {

          rm(p, x)

          gc()

        }

        readr::write_lines(paste0("predict finished in "
                                  , round(difftime(Sys.time(), pred_start, units = "mins"), 2)
                                  , " minutes"
                                  )
                           , file = log_file
                           , append = TRUE
                           )

      }

      run <- if(file.exists(thresh_file)) force_new else TRUE

      if(all(apply_thresh, run)) {

        ## thresh ------
        thresh_start <- Sys.time()

        m <- paste0("thresh for ", this_taxa)

        message(m)

        readr::write_lines(m
                           , file = log_file
                           , append = TRUE
                           )

        thresh <- mod$e[[1]]@thresholds$max_spec_sens

        terra::app(terra::rast(pred_file)
                   , \(i) i > thresh
                   , filename = thresh_file
                   , overwrite = TRUE
                   , wopt = list(datatype = "INT1U"
                                 , names = this_taxa
                                 )
                   )

        readr::write_lines(paste0("thresh done in "
                                  , round(difftime(Sys.time(), thresh_start, units = "mins"), 2)
                                  , " minutes"
                                  )
                           , file = log_file
                           , append = TRUE
                           )

        if(do_gc) {

          gc()

        }

      }

      # finish ----
      readr::write_lines(paste0("predict finished. elapsed time: "
                                , round(difftime(Sys.time(), start_time, units = "mins"), 2)
                                , " minutes"
                                )
                         , file = log_file
                         , append = TRUE
                         )

    } else {

      m <- paste0("predict did not run as:\n "
                  , if(prep$abandoned) "prep abandoned\n "
                  , if(! prep$finished) "prep did not finish\n "
                  , if(any(prep$abandoned, ! prep$finished)) paste0("prep log:\n  ", paste0(prep$log, collapse = "\n  "), "\n ")
                  , if (! full_run$finished) paste0("full run did not finish\n full run log:\n  ", paste0(full_run$log, collapse = "\n  "), "\n ")
                  , if(file.exists(thresh_file)) paste0("threshold file: ", thresh_file, " already exists and force_new was not TRUE")
                  , "\n"
                  )

      readr::write_lines(m
                         , file = log_file
                         , append = TRUE
                         )

    }

    if(do_gc) {

      stuff <- ls()

      delete_stuff <- stuff[!grepl("_file$", stuff)]

      rm(list = delete_stuff)

      gc()

    }

    res <- c(if(file.exists(pred_file)) list(pred = pred_file)
             , if(file.exists(thresh_file)) list(thresh = thresh_file)
             )

    return(res)

  }

