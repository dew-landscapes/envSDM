
#' Predict from  SDM
#'
#' @param this_taxa Character. Name of taxa. Used to name outputs. If `NULL`,
#' this will be `basename(dirname(out_dir))`.
#' @param prep_dir Character. Name of directory containing: `prep.rds` (created
#' with `envSDM::prep_sdm()`)
#' @param tune_dir Character. Name of directory containing `tune.rds`, created
#' with `envSDM::tune_sdm()`. Note that any `tune.rds` can be used but only the
#' model in the first row will be used, thus more usually this `tune.rds` will
#' have been created directly by `envSDM::run_full_sdm()`
#' @param out_dir Character. Name of directory into which `.tif`s will be saved.
#' Will be created if it does not exist.
#' @param predictors Character. Vector of paths to predictor `.tif` files.
#' @param is_env_pred Logical. Does the naming of the directory and files in
#' `predictors` follow the pattern required by `envRaster::parse_env_tif()`?
#' @param terra_options Passed to `terra::terraOptions()`. e.g. list(memfrac = 0.6)
#' @param doClamp Passed to `terra::predict()` (which then passes as `...` to
#' `fun`). Possibly orphaned from older envSDM?
#' @param limit_to_mcp Logical. If `predict_boundary` exists within `prep` and
#' `limit_to_mcp == TRUE`, an output raster (`mask.tif`) will be created within
#' `predict_boundary` using `terra::mask()`. Irrespective of `limit_to_mcp`,
#' `full.tif` is always created at the full extent of the predictors. Thus all
#' `mask.tif` files can be 'stacked' as they have the same extent. If needed,
#' limiting the predictions for a taxa to its predict boundary can then be done
#' via `terra::trim(mask.tif)`.
#' @param apply_thresh Logical. If `TRUE`, an output raster `thresh.tif` will be
#' created using the threshold `mod$e[[1]]@thresholds$max_spec_sens`
#' @param force_new Logical. If outputs already exist, should they be remade?
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful to keep RAM use down when running SDMs for many, many taxa,
#' especially if done in parallel.
#' @param check_tifs Logical. Check if any output `.tif` files error on
#' `terra::rast()` and delete them if they do. Useful after a crash during
#' predict.
#' @param ... Passed to terra::predict. e.g. use for wopt = list()
#'
#' @return `invisible(NULL)`. Output .tif, .log, and optional .png, written to
#' `out_dir`
#' @export
#'
#' @example inst/examples/predict_sdm_ex.R
#'
  predict_sdm <- function(this_taxa = NULL
                          , prep_dir
                          , tune_dir = NULL
                          , out_dir = NULL
                          , predictors = NULL
                          , is_env_pred = TRUE
                          , terra_options = NULL
                          , doClamp = TRUE
                          , limit_to_mcp = TRUE
                          , apply_thresh = TRUE
                          , force_new = FALSE
                          , do_gc = FALSE
                          , check_tifs = TRUE
                          , ...
                          ) {

    this_taxa <- basename(prep_dir)

    if(is.null(tune_dir)) tune_dir <- prep_dir
    if(is.null(out_dir)) out_dir <- tune_dir


    # files -----
    ## existing
    prep_file <- fs::path(prep_dir, "prep.rds")
    prep_log <- fs::path(prep_dir, "prep.log")
    tune_file <- fs::path(tune_dir, "tune.rds")

    ## new
    pred_file <- fs::path(out_dir, "full.tif")
    mask_file <- fs::path(out_dir, "mask.tif")
    thresh_file <- fs::path(out_dir, "thresh.tif")
    pred_log <- fs::path(out_dir, "pred.log")


    # test tifs ------
    # need to test before running, in case the test deletes an incomplete .tif
    if(all(check_tifs, dir.exists(out_dir))) {

      safe_rast <- purrr::safely(terra::rast)

      test_files <- fs::dir_ls(out_dir
                               , regexp = "tif$|nc$"
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
    log_text <- paste0(readLines(prep_log), collapse = " ")

    run <- all(file.exists(prep_file) # need a prep file
               , !grepl("SDM abandoned", log_text) # enough records
               , grepl("prep end", log_text) # prep finished
               , if(file.exists(thresh_file)) force_new else TRUE # output not already exists (unless force_new)
               , file.exists(tune_file) # need a tune to predict from
               )

    if(run) {

      fs::dir_create(out_dir)

      prep <- rio::import(prep_file)

      # start timer ------
      pred_timer <- envFunc::timer(process = "predict start"
                                   , file = "predict"
                                   , time_df = NULL
                                   , log = pred_log
                                   , write_log = TRUE
                                   )

      mod <- rio::import(fs::path(tune_dir, "tune.rds")
                         , setclass = "tibble"
                         )

      if(nrow(mod) > 1) {

        warning(nrow(mod)
                , " tunes present. Using the tune from the first row."
                )

      }

      ## limit -----
      if(!exists("predict_boundary", where = prep)) {

        prep$pred_limit <- FALSE

      } else {

        prep$pred_limit <- limit_to_mcp

      }

      pred_limit <- prep$pred_limit

      # pred -------

      run <- if(file.exists(pred_file)) force_new else TRUE

      if(all(!is.null(predictors), run)) {

        gc()

        algo <- mod$algo[[1]]

        if(algo == "maxnet") requireNamespace("maxnet", quietly = TRUE)
        if(algo == "rf") requireNamespace("randomForest", quietly = TRUE)
        if(algo == "env") requireNamespace("predicts", quietly = TRUE)

        ## predictors -----
        if(is_env_pred) {

          pred_names <- envRaster::name_env_tif(tibble::tibble(path = predictors), parse = TRUE) %>%
            dplyr::mutate(name = paste0(season,"__", layer)) %>%
            dplyr::pull(name)

          x <- terra::rast(predictors)
          names(x) <- pred_names

        } else {

          x <- terra::rast(predictors)

        }

        safe_predict <- purrr::safely(terra::predict)

        if(!is.null(terra_options)) {

          do.call(terra::terraOptions
                  , args = terra_options
                  )

        }


        # full -----
        message("full predict for ", this_taxa)

        p <- safe_predict(object = x
                          , model = mod$m[[1]]
                          , type = if(algo == "rf") "prob" else "cloglog"
                          , clamp = doClamp
                          , na.rm = TRUE
                          , index = if(algo == "rf") "1" else NULL
                          , names = this_taxa
                          , filename = pred_file
                          , overwrite = TRUE
                          , ...
                          )

        if(!is.null(p$error)) {

          pred_timer <- envFunc::timer("full"
                                       , notes = as.character(p$error)
                                       , time_df = pred_timer
                                       , write_log = TRUE
                                       )

        }

        if(do_gc) {

          rm(p, x)

          gc()

        }

        pred_timer <- envFunc::timer("full"
                            , time_df = pred_timer
                            )

      }

      run <- if(file.exists(mask_file)) force_new else TRUE

      if(all(pred_limit, run, file.exists(pred_file))) {

        # mask ------
        message("mask for ", this_taxa)

        terra::mask(terra::rast(pred_file)
                    , terra::vect(prep$predict_boundary)
                    , filename = mask_file
                    , overwrite = TRUE
                    , wopt = list(datatype = "INT2S"
                                  , scale = gdalcubes::pack_minmax(min = 0, max = 1)$scale
                                  , offset = gdalcubes::pack_minmax(min = 0, max = 1)$offset
                                  , NAflag = gdalcubes::pack_minmax(min = 0, max = 1)$nodata
                                  )
                    )

        pred_timer <- envFunc::timer("mask"
                            , time_df = pred_timer
                            )

        if(do_gc) {

          gc()

        }

      }

      # thresh------

      run <- if(file.exists(thresh_file)) force_new else TRUE

      if(all(apply_thresh, run, file.exists(pred_file), file.exists(tune_file))) {

        do_thresh <- if(file.exists(thresh_file)) force_new else TRUE

        if(do_thresh) {

          thresh <- mod$e[[1]]@thresholds$max_spec_sens

          file_to_thresh <- if(pred_limit) {

            mask_file

          } else {

            pred_file

          }

          message("apply threshold for ", this_taxa)

          terra::app(terra::rast(file_to_thresh)
                     , \(i) i > thresh
                     , filename = thresh_file
                     , overwrite = TRUE
                     , wopt = list(datatype = "INT1U")
                     )

          pred_timer <- envFunc::timer("threshold"
                              , time_df = pred_timer
                              )

          if(do_gc) {

            gc()

          }

        }

      }

      pred_timer <- envFunc::timer(process = "predict end"
                          , time_df = pred_timer
                          )

    }

    return(invisible(NULL))

  }

