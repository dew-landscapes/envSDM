
#' Run an SDM using no cross validation and previously established tune arguments
#'
#' @param this_taxa Character. Name of taxa. Used to name outputs. If `NULL`,
#' this will be `basename(dirname(out_dir))`.
#' @param out_dir Character. Name of directory containing previous tunes and
#' into which results will be saved.
#' @param use_metric Character. Which metric to use to find the 'best' tune
#' arguments from previous tuning results? Default is `combo`, the product of
#' `auc_po`, `CBI_rescale` and `IMAE`. `use_metric` must be `combo` or have been
#'  used in the use_metrics argument to `tune_sdm()`.
#' @param force_new Logical. If outputs already exist, should they be remade?
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param save_to Character. Name of path to save results. Defaults to
#' `fs::path(out_dir, use_metric)`
#' @param ... Passed to `tune_sdm()`
#'
#' @return Character path to output .rds file. `tune.rds` saved into `save_to`
#' directory. log written. `tune.rds` is a data frame but performs poorly, due
#' to list columns, if not imported as a tibble (e.g. via
#' `rio::import("tune.rds", setclass = "tibble"))`)
#' @export
#'
#' @examples inst/examples/predict_sdm_ex.R
  run_full_sdm <- function(this_taxa = NULL
                           , out_dir
                           , use_metric = "combo"
                           , force_new = FALSE
                           , save_to = fs::path(out_dir, use_metric)
                           , do_gc = FALSE
                           , ...
                           ) {

    message(paste0("run full model for ", basename(out_dir)))

    # prep -----

    ## files ------
    # existing
    prep_file <- fs::path(out_dir, "prep.rds")

    prep_log <- fs::path(out_dir, "prep.log")

    ## new
    tune_file <- fs::path(save_to, "tune.rds")
    full_sdm_log <- fs::path(save_to, "full_sdm.log")

    # run?-----
    prep_log_text <- paste0(readLines(prep_log), collapse = " ")

    run <- all(file.exists(prep_file)
               , !grepl("SDM abandoned", prep_log_text)
               , grepl("prep end", prep_log_text)
               )

    if(run) {

      run_tune <- if(file.exists(tune_file)) force_new else TRUE

      if(run_tune) {

        # to create
        fs::dir_create(save_to)

        # start timer ------
        full_sdm_timer <- envFunc::timer(process = "full SDM start"
                                         , file = "full"
                                         , time_df = NULL
                                         , log = full_sdm_log
                                         , write_log = TRUE
                                         )

        # tune_args------
        tune_args <- rio::import(fs::path(out_dir, "evaluation.csv")
                                 , setclass = "tibble"
                                 ) %>%
          dplyr::mutate(filter_col = !!rlang::ensym(use_metric)) %>%
          dplyr::filter(filter_col == max(filter_col))

        if(length(tune_args)) {

          # mod --------

          tune_sdm(out_dir = out_dir
                   , algo = tune_args$algo
                   , fc = tune_args$fc
                   , rm = tune_args$rm
                   , trees = tune_args$trees
                   , mtry = tune_args$mtry
                   , nodesize = tune_args$nodesize
                   , keep_model = TRUE
                   , best_run = TRUE
                   , do_gc = do_gc
                   , save_to = save_to
                   , ...
                   )

          full_sdm_timer <- envFunc::timer("full SDM end"
                                           , time_df = full_sdm_timer
                                           )

          if(do_gc) {

            stuff <- ls()

            delete_stuff <- stuff[stuff != "tune_file"]

            rm(list = delete_stuff)

            gc()

          }

        } else {

          full_sdm_timer <- envFunc::timer("warning"
                                           , notes = "No tune results found"
                                           , time_df = full_sdm_timer
                                           )

          full_sdm_timer <- envFunc::timer("full sdm end"
                                           , timer_df = full_sdm_timer
                                           )

        }

      }

    }

    res <- if(file.exists(tune_file)) tune_file else NULL

    return(res)

  }

