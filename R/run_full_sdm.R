
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
#' a path to the saved file. If `out_dir` is a valid path, the 'full
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
    # everything should go through here, at least as far as blah_file.rds
    if(is.character(out_dir)) {

      fs::dir_create(out_dir)

      if(dir.exists(out_dir)) {

        full_run_file <- fs::path(out_dir
                                  , "full_run.rds"
                                  )

        if(all(file.exists(full_run_file), ! force_new)) {

          safe_import <- purrr::safely(rio::import)

          full_run <- safe_import(full_run_file, trust = TRUE)

          if(is.null(full_run$error)) full_run <- full_run$result else {

            # remove the full_run file if it can't be opened
            fs::file_delete(full_run_file)

            rm(full_run)

          }

        }

      } else stop("can't create out_dir")

    }

    ## prep -------
    if(! "list" %in% class(prep)) prep <- rio::import(prep, trust = TRUE)

    ## tune ---------
    if(! "list" %in% class(tune)) tune <- rio::import(tune, trust = TRUE)

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
      row_id <- 1

      tune_args <- tune$tune_mean %>%
        dplyr::filter(tunes == max(tunes)) |>
        dplyr::arrange(desc(!!rlang::ensym(use_metric)))


      full_run_tune <- NULL

      if(nrow(tune_args)) {

        # mod --------

        while(is.null(full_run_tune$tune_mean)) {

          full_run_tune <- tune_sdm(prep = prep
                                    , out_dir = out_dir
                                    , return_val = "object"
                                    , algo = tune_args$algo[[row_id]]
                                    , fc = tune_args$fc[[row_id]]
                                    , rm = tune_args$rm[[row_id]]
                                    , trees = tune_args$trees[[row_id]]
                                    , mtry = tune_args$mtry[[row_id]]
                                    , nodesize = tune_args$nodesize[[row_id]]
                                    , keep_model = TRUE
                                    , best_run = TRUE
                                    , do_gc = do_gc
                                    , force_new = TRUE
                                    , ...
                                    )

          if(is.null(full_run_tune$tune_mean)) row_id = row_id + 1

        }

        if(row_id > 2) {

          message("error using tune arguments: "
                  , tune_args$tune_args[row_id - 1]
                  , ". trying arguments with next best "
                  , use_metric
                  , " value"
                  )

        }

        full_run$finished <- TRUE

        full_run <- c(full_run, full_run_tune)

      }

    }

    # save -------
    # export before gc()
    rio::export(full_run, full_run_file)

    if(do_gc) {

      stuff <- ls()

      delete_stuff <- stuff[! stuff %in% c(return_val, "return_val")]

      rm(list = delete_stuff)

      gc()

    }

    res <- if(return_val == "full_run") full_run else full_run_file

    return(res)

  }

