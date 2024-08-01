
#' Predict an SDM
#'
#' @param this_taxa Character. Name of taxa. Used to name outputs. If `NULL`,
#' this will be `basename(dirname(out_dir))`.
#' @param out_dir Character. Name of directory containing model to predict from
#' and into which results will be saved.
#' @param predictors Character. Directory in which predictor `.tif` files live.
#' `envRaster::name_env_tif()` is used internally to parse so should follow that
#' naming convention. :( (might work otherwise, but seems unlikely).
#' @param terra_options Passed to `terra::terraOptions()`. e.g. list(memfrac = 0.6)
#' @param doClamp Passed to `terra::predict()` (which then passes as `...` to
#' `fun`). Possibly orphaned from older envSDM?
#' @param limit_to_mcp Logical. If `predict_boundary` exists within `prep` and
#' `limit_to_mcp == TRUE`, the output raster will be limted to within
#' `predict_boundary`.
#' @param apply_thresh Logical. If `TRUE`, an output raster `thresh.tif` will be
#' created using the threshold `mod$e[[1]]@thresholds$max_spec_sens`
#' @param force_new Logical. If outputs already exist, should they be remade?
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param check_tifs Logical. Check that any output `.tif` files error on
#' `terra::rast()` and delete them if they do. Useful after crash during predict.
#' @param ... Not used.
#'
#' @return `invisible(NULL)`. Output `.tif` files are created.
#' @export
#'
#' @examples
#'
#'
  predict_sdm <- function(this_taxa = NULL
                          , out_dir
                          , predictors = NULL
                          , terra_options = NULL # list(memfrac = 0.6)
                          , doClamp = TRUE
                          , limit_to_mcp = TRUE
                          , apply_thresh = TRUE
                          , force_new = FALSE
                          , do_gc = FALSE
                          , check_tifs = TRUE
                          , ...
                          ) {

    this_taxa <- basename(dirname(out_dir))

    # files -----
    ## existing
    prep_file <- fs::path(dirname(out_dir), "prep.rds")
    tune_file <- fs::path(out_dir, "tune.rds")
    prep_log <- fs::path(dirname(out_dir), "prep.log")

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

      prep <- rio::import(prep_file)

      # start timer ------
      pred_timer <- envFunc::timer(process = "predict start"
                          , file = "predict"
                          , time_df = NULL
                          , log = pred_log
                          , write_log = TRUE
                          )

      mod <- rio::import(fs::path(out_dir, "tune.rds")
                         , setclass = "tibble"
                         )

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
        pred_names <- envRaster::name_env_tif(tibble::tibble(path = predictors), parse = TRUE) %>%
          dplyr::mutate(name = paste0(season,"__", layer)) %>%
          dplyr::pull(name)

        x <- terra::rast(predictors)
        names(x) <- pred_names

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
                          , wopt = list(datatype = "INT2S"
                                        , scale = 0.00001525925
                                        , offset = 0.5
                                        , NAflag = -32768
                                        )
                          )

        if(!is.null(p$error)) {

          readr::write_lines(p$error
                             , file = gsub("tif$", "log", pred_file)
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
                                  , scale = 0.00001525925
                                  , offset = 0.5
                                  , NAflag = -32768
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

