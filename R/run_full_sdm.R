
#' Run an SDM using no cross validation and previously established tune arguments
#'
#' @param prep Character or named list. If character, the path to an existing
#' `prep.rds`. Otherwise, the result of a call to prep_sdm with return_val =
#' "object"
#' @param tune Character or named list. If character, the path to an existing
#' `tune.rds`. Otherwise, the result of a call to tune_sdm with return_val =
#' "object"
#' @param out_dir FALSE or character. If FALSE the result of `run_full_sdm()`
#' will be saved to a temporary folder. If character, a file 'tune.rds' will be
#' created at the path defined by `out_dir`.
#' @param return_val Character: "object" or "path". Both return a named list. In
#' the case of "path" the named list is simply list(full_run = out_dir). Will
#' be set to "object" if `out_dir` is FALSE.
#' @param use_metric Character. Which metric to use to find the 'best' tune
#' arguments from previous tuning results? Default is `combo`, the product of
#' `auc_po`, `CBI_rescale` and `IMAE`. `use_metric` must be `combo` or have been
#'  used in the use_metrics argument to `tune_sdm()`.
#' @param force_new Logical. If outputs already exist, should they be remade?
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param ... Passed to `tune_sdm()`
#'
#' @return If `return_val` is "object" a named list. If `return_val` is "path"
#' a named list `list(prep = out_dir)`. If `out_dir` is a valid path, the 'full
#' result' (irrespective of `return_val`) is also saved to
#' `fs::path(out_dir, "prep.rds")`. The 'full result' is a named list with
#' elements:
#'
#' @export
#'
#' @examples inst/examples/predict_sdm_ex.R
  run_full_sdm <- function(prep
                           , tune
                           , out_dir
                           , return_val = "path"
                           , use_metric = "combo"
                           , force_new = FALSE
                           , do_gc = FALSE
                           , ...
                           ) {

    # setup -------
    ## return ------
    return_val <- if(any(isFALSE(out_dir), return_val == "object")) "full_run" else "full_run_file"

    if(isFALSE(out_dir)) out_dir <- tempfile()

    ## out_dir ------
    if(is.character(out_dir)) {

      fs::dir_create(out_dir)

      if(dir.exists(out_dir)) {

        full_run_file <- fs::path(out_dir
                              , "full_run.rds"
                              )

        if(file.exists(full_run_file)) {

          full_run <- rio::import(full_run_file)

        }

      } else stop("can't create out_dir")

    }

    ## prep -------
    if(! "list" %in% class(prep)) prep <- rio::import(prep)

    ## tune ---------
    if(! "list" %in% class(tune)) tune <- rio::import(tune)

    ## full_run -------
    if(!exists("full_run", inherits = FALSE)) full_run <- list(finished = FALSE)

    # run?-----
    run <- all(!prep$abandoned
               , prep$finished
               , tune$finished
               , if(full_run$finished) force_new else TRUE
               )

    if(run) {

      # tune_args------
      tune_args <- tune$tune_mean %>%
        dplyr::filter(!!rlang::ensym(use_metric) == max(!!rlang::ensym(use_metric), na.rm = TRUE))

      if(length(tune_args)) {

        # mod --------

        full_run_tune <- tune_sdm(prep = prep
                                  , out_dir = out_dir
                                  , return_val = "object"
                                  , algo = tune_args$algo
                                  , fc = tune_args$fc
                                  , rm = tune_args$rm
                                  , trees = tune_args$trees
                                  , mtry = tune_args$mtry
                                  , nodesize = tune_args$nodesize
                                  , keep_model = TRUE
                                  , best_run = TRUE
                                  , do_gc = do_gc
                                  , ...
                                  )

        full_run <- c(full_run, full_run_tune)

      }

    }

    # save -------
    # export before gc()
    full_run$finished <- TRUE
    rio::export(full_run, full_run_file)

    if(do_gc) {

      stuff <- ls()

      delete_stuff <- stuff[! stuff %in% c(return_val, "return_val")]

      rm(list = delete_stuff)

      gc()

    }

    res <- if(return_val == "full_run") get("full_run") else list(full_run_file = get("full_run_file"))

    return(res)

  }

