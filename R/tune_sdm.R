
#' Tune, and evaluate, species distribution models
#'
#' @param prep Character or named list. If character, the path to an existing
#' `prep.rds`. Otherwise, the result of a call to prep_sdm with return_val =
#' "object"
#' @param out_dir FALSE or character. If FALSE the result of tune_sdm will be
#' saved to a temporary folder. If character, a file 'tune.rds' will be created
#' at the path defined by out_dir.
#' @param return_val Character: "object" or "path". Both return a named list. In
#' the case of "path" the named list is simply list(tune = out_dir). Will be set
#' to "object" if `out_dir` is FALSE.
#' @param algo Character. Name of algorithm to use.
#' @param fc Character. Used to generate levels of `classes` argument to
#' `maxnet::maxnet()` that are tuned.
#' @param limit_p `TRUE`, `FALSE` or number of predictor variables above which
#' to limit the use of `p` in the classes argument used in `maxnet::maxnet()`.
#' Useful with many predictor variables when it becomes unwieldy to generate
#' interactions for all predictors.
#' @param rm Numeric. Used to generate levels of `regmult` argument to
#' `maxnet::maxnet()` that are tuned.
#' @param trees Used to generate the levels of `ntree` argument to
#' `randomForest::randomForest()` that are tuned. `TRUE` (tune with default
#' `trees`), `FALSE` (don't tune `trees`) or numeric (the `trees` values to tune
#'  with).
#' @param mtry Used to generate the levels of `mtry` argument to
#' `randomForest::randomForest()` that are tuned. `TRUE` (tune with sensible guesses for
#' `mtry`), `FALSE` (only use default `randomForest::randomForest()` `mtry`) or
#' numeric (the `mtry` values to tune with).
#' @param limit_spat_mtry Numeric. If `mtry` is `TRUE` and if using spatial
#' cross validation, the values of `mtry` to tune will be limited to less than
#' or equal to `limit_spat_mtry`.
#' @param nodesize Used to generate the levels of `nodesize` argument to
#' `randomForest::randomForest()` that are tuned. `TRUE` (tune with default
#' `nodesize`), `FALSE` (only use default `randomForest::randomForest()`
#' `nodesize`) or numeric (the `nodesize` values to tune with).
#' @param keep_model Logical. If `TRUE` the model results will be appended as a
#' list column in the returned tibble (as column `m`)
#' @param best_run Logical. If `TRUE` this alters the behaviour of the
#' `tune_sdm()` by, well, not tuning. :). Sets all blocks to the same value so
#' no cross-validation.
#' @param metrics_df Dataframe. Defines which metrics to use when deciding on
#' 'good' SDMs.
#' @param use_metrics Character. Vector of values in metrics_df$metric to use
#' when finding the 'best' model.
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param force_new Logical. If outputs already exist, should they be remade?
#' @param ... Passed to `evaluate_sdm()`. e.g. thresholds for use in
#' `predicts::pa_evaluate()` (as `tr` argument, although if used, the values of
#' the `thresholds` element of the `pa_ModelEvaluation` object returned by
#' `predicts::pa_evaluate()` will be limited to the values in `tr`).
#'
#' @return If `return_val` is "object" a named list. If `return_val` is "path"
#' a named list `list(prep = out_dir)`. If `out_dir` is a valid path, the 'full
#' result' (irrespective of `return_val`) is also saved to
#' `fs::path(out_dir, "prep.rds")`. The 'full result' is a named list with
#' elements:
#'
#' @export
#'
#' @example inst/examples/tune_sdm_ex.R
  tune_sdm <- function(prep
                       , out_dir = FALSE
                       , return_val = "path"
                       , algo = c("all", "maxnet", "bioclim", "envelope", "rf") # envelope = bioclim
                       , fc = "auto_feature" # maxnet tune
                       , limit_p = FALSE # T, F or number of preds above which to limit p
                       , rm = seq(1, 6, 0.5) # maxnet tune
                       , trees = c(500, 1000, 2000) # T, F or numeric
                       , mtry = TRUE # T, F or numeric
                       , limit_spat_mtry = 4
                       , nodesize = c(1, 2) # T, F or numeric
                       , keep_model = FALSE
                       , best_run = FALSE
                       , metrics_df = envSDM::sdm_metrics
                       , use_metrics = c("auc_po", "CBI_rescale", "IMAE")
                       , do_gc = TRUE
                       , force_new = FALSE
                       , ...
                       ) {

    # setup -------
    ## return ------
    return_val <- if(any(isFALSE(out_dir), return_val == "object")) "tune" else "tune_file"

    if(isFALSE(out_dir)) out_dir <- tempfile()

    ## log file ------
    log_file <- fs::path(out_dir
                         , if(!best_run) "tune.log" else "full_run.log"
                         )

    ## out_dir ------
    if(is.character(out_dir)) {

      fs::dir_create(out_dir)

      if(dir.exists(out_dir)) {

        tune_file <- fs::path(out_dir
                              , if(!best_run) "tune.rds" else "full_run.rds"
                              )

        if(file.exists(tune_file)) {

          tune <- rio::import(tune_file, trust = TRUE)

        }

      } else stop("can't create out_dir")

    }

    ## prep -------
    if(! "list" %in% class(prep)) prep <- rio::import(prep, trust = TRUE)

    ## tune ---------
    if(!exists("tune", inherits = FALSE)) tune <- list(finished = FALSE)

    # run?-----
    run <- all(!prep$abandoned
               , prep$finished
               , if(tune$finished) force_new else TRUE
               )

    if(run) {

      this_taxa <- prep$this_taxa

      message("tuning "
              , this_taxa
              , "\nout_dir is "
              , out_dir
              )

      # tune -----
      if(exists("blocks", where = prep)) {

        # start timer ------
        start_time <- Sys.time()

        readr::write_lines(paste0("\n\n"
                                  , this_taxa
                                  , "\nbest_run = "
                                  , best_run
                                  , "\ntune start "
                                  , start_time
                                  )
                           , file = log_file
                           , append = TRUE
                           )


          ## setup ------

          nobs <- nrow(prep$blocks[prep$blocks$pa == 1,])

          ## start data frame -----

          if(best_run) {

            old_prep_block <- prep$blocks$block

            prep$blocks$block <- 1

          }

          single_block <- if(length(unique(prep$blocks[prep$blocks$pa == 1,]$block)) == 1) TRUE else FALSE

          preds <- prep$reduce_env$env_cols[prep$reduce_env$env_cols %in% prep$reduce_env$keep]

          readr::write_lines(paste0(length(preds)
                                    , " out of "
                                    , length(prep$reduce_env$env_cols)
                                    , " variables will be used"
                                    )
                             , file = log_file
                             , append = TRUE
                             )

          p <- prep$blocks$pa

          data <- prep$blocks %>%
            dplyr::select(tidyselect::any_of(preds))

          # start_df --------
          start_df <- tibble::tibble(k = sort(unique(prep$blocks$block))) %>%
            dplyr::mutate(ids = if(any(best_run, single_block)) purrr::map(k, ~prep$blocks$block == .) else purrr::map(k, ~ prep$blocks$block != .)
                          , pa_train = purrr::map(ids, \(a) p[a])
                          , data_train = purrr::map(ids, \(a) data[a,])
                          , p_data_test = purrr::map(k
                                                     , \(a) prep$blocks %>%
                                                       dplyr::filter(p == 1
                                                                     , block == a
                                                                     ) %>%
                                                       dplyr::select(tidyselect::any_of(preds)) %>%
                                                       as.data.frame()
                                                     )
                          , a_data_test = purrr::map(k
                                                     , \(a) prep$blocks %>%
                                                       dplyr::filter(p == 0
                                                                     , block == a
                                                                     ) %>%
                                                       dplyr::select(tidyselect::any_of(preds)) %>%
                                                       as.data.frame()
                                                     )
                          , n_p_test = purrr::map_dbl(p_data_test
                                                      , nrow
                                                      )
                          , n_p_train = purrr::map_dbl(pa_train, \(x) sum(x == 1))
                          , n_a_train = purrr::map2_dbl(data_train
                                                        , n_p_train
                                                        , \(x, y) nrow(x) - y
                                                        )
                          # 'hack' to ensure min of p's and a's is used (as n_p_train)
                          , n_p_train = purrr::map2_dbl(n_p_train
                                                        , n_a_train
                                                        , \(x, y) min(x, y)
                                                        )
                          ) %>%
            dplyr::filter(n_p_test > 0)

          if(best_run) {

            # testing set --------
            start_df <- start_df %>%
              dplyr::mutate(p_data_test = list(prep$testing %>%
                                                 dplyr::filter(pa == 1) %>%
                                                 dplyr::select(tidyselect::any_of(preds)) %>%
                                                 as.data.frame()
                                               )
                            , a_data_test = list(prep$testing %>%
                                                   dplyr::filter(pa == 0) %>%
                                                   dplyr::select(tidyselect::any_of(preds)) %>%
                                                   as.data.frame()
                                                 )
                            , n_a_train = sum(prep$testing$pa == 1)
                            )


          }

          ## tune maxnet------

          if(any(c("all", "maxnet") %in% algo)) {

            run <- if(exists("tune_maxnet", tune, inherits = FALSE)) force_new else TRUE

            if(run) {

              message("maxnet tune")
              start_maxnet <- Sys.time()

              ### fc --------
              # check auto_features
              if (any(c(fc == "auto_feature", c("Q", "P", "T", "H") %in% fc))) {

                if(any(grepl("auto_feature", fc))) fc <- c("L", "Q", "H", "LQ", "QH", "LQP", "QHP", "LQHP")

                # set fc based on number of obs
                if (nobs <= 10) {
                  fc <- fc[fc %in% "L"] # only linear if less than 10 recs
                } else if (nobs <= 15) {
                  fc <- fc[fc %in% c("L", "Q", "LQ")]
                } else if (nobs < 80) {
                  fc <- fc[fc %in% c("L", "Q", "LQ", "LH", "QH", "LQH")]
                }
              }

              # make sure any fc's requested are unique and H and T aren't both specified
              fc <- unique(fc)

              # if threshold AND hinge requested, keep hinge only
              if (all(any(grepl(pattern = "H", fc)), any(grepl(pattern = "T", fc)))) {
                # remove T
                fc <- unique(gsub("T", "", fc))
                # remove any empty
                fc <- fc[nchar(fc) > 0]
              }

              ## limit_p ------
              # if limit_p, ensure p not in fc
              if(!isFALSE(limit_p)) {

                # Only limit if too many preds, or limit_p is TRUE
                if(is.numeric(limit_p)) {

                  remove_p <- length(preds) > limit_p

                } else if(isTRUE(limit_p)) {

                  remove_p <- TRUE

                } else {

                  stop("limit_p must be TRUE, FALSE or numeric. current vaule is: "
                       , limit_p
                       )

                }

                if(remove_p) {
                  # remove P
                  fc <- unique(gsub("P", "", fc))
                  # remove any empty
                  fc <- fc[nchar(fc) > 0]
                }

              }

              safe_maxnet <- purrr::safely(maxnet::maxnet)

              tune_maxnet <- start_df %>%
                dplyr::cross_join(tibble::tibble(fc = tolower(fc))) %>%
                dplyr::cross_join(tibble::tibble(rm = rm)) %>%
                dplyr::mutate(tune_args = paste0("fc: ", fc
                                                 , ". rm: ", rm
                                                 )
                              , m = purrr::pmap(list(pa_train
                                                     , data_train
                                                     , fc
                                                     , rm
                                                     )
                                              , \(a, b, c, d) safe_maxnet(p = a
                                                                          , data = b
                                                                          , f = maxnet::maxnet.formula(p = a
                                                                                                        , data = b
                                                                                                        , classes = c
                                                                                                       )
                                                                          , regmult = d
                                                                          )
                                              )
                              )

              probs <- tune_maxnet %>%
                dplyr::mutate(err = purrr::map(m, "error")) %>%
                dplyr::filter(purrr::map_lgl(err, \(x) !is.null(x)))

              if(nrow(probs)) {

                purrr::map2(probs$err
                            , probs$tune_args
                            , \(a, b) {
                              readr::write_lines(paste0("error with tune_args "
                                                        , b
                                                        , ": "
                                                        , as.character(a)
                                                        )
                                                 , file = log_file
                                                 , append = TRUE
                                                 )
                              }
                           )

              }

              tune$tune_maxnet <- tune_maxnet %>%
                dplyr::mutate(m = purrr::map(m, "result")) %>%
                dplyr::filter(purrr::map_lgl(m, \(x) !is.null(x))) %>%
                dplyr::mutate(e = purrr::pmap(list(m
                                                   , p_data_test
                                                   , a_data_test
                                                   )
                                              , \(a, b, c) evaluate_sdm(a, b, c
                                                                        , ...
                                                                        , do_gc = do_gc
                                                                        )
                                              )
                              )

              null_e <- purrr::map_lgl(tune$tune_maxnet$e, \(x) is.null(x))

              if(sum(null_e)) {

                null_e_tune <- envFunc::vec_to_sentence(unique(tune$tune_maxnet$tune_args[null_e]))
                blocks <- envFunc::vec_to_sentence(unique(tune$tune_maxnet$k[null_e]))

                readr::write_lines(paste0("warning. tune_args: "
                                          , null_e_tune
                                          , " failed to predict ("
                                          , sum(null_e)
                                          , " out of "
                                          , nrow(tune$tune_maxnet)
                                          , " total tunes) in block(s) "
                                          , blocks
                                          )
                                 , file = log_file
                                 , append = TRUE
                                 )

              }

              tune$tune_maxnet <- tune$tune_maxnet %>%
                dplyr::filter(purrr::map_lgl(e, \(x) ! is.null(x))) %>%
                {if(keep_model) (.) %>% dplyr::select(! dplyr::where(is.list), m, e) else (.) %>% dplyr::select(! dplyr::where(is.list), e)}

              readr::write_lines(paste0("maxnet tune finished in: "
                                        , round(difftime(Sys.time(), start_maxnet, units = "mins"), 2)
                                        , " minutes"
                                        )
                                 , file = log_file
                                 , append = TRUE
                                 )

            }

          }


          ## tune envelope -------

          if(any(c("all", "envelope", "bioclim", "env") %in% algo)) {

            run <- if(exists("tune_envelope", tune, inherits = FALSE)) force_new else TRUE

            if(run) {

              message("envelope tune")

              start_envelope <- Sys.time()

              tune$tune_envelope <- start_df %>%
                dplyr::mutate(m = purrr::map(data_train
                                             , function(x) predicts::envelope(x = x)
                                             )
                              , tune_args = "none"
                              , e = purrr::pmap(list(m
                                                     , p_data_test
                                                     , a_data_test
                                                     )
                                                , \(a, b, c) evaluate_sdm(a, b, c
                                                                          , ...
                                                                          , do_gc = do_gc
                                                                          )
                                                )
                              ) %>%
                {if(keep_model) (.) %>% dplyr::select(! dplyr::where(is.list), m, e) else (.) %>% dplyr::select(! dplyr::where(is.list), e)}

              readr::write_lines(paste0("envelope tune finished in: "
                                        , round(difftime(Sys.time(), start_envelope, units = "mins"), 2)
                                        , " minutes"
                                        )
                                 , file = log_file
                                 , append = TRUE
                                 )

            }

          }

          ## tune rf --------

          if(any(c("all", "rf", "randomForest") %in% algo)) {

            run <- if(exists("tune_rf", tune, inherits = FALSE)) force_new else TRUE

            if(run) {

              message("rf tune")

              start_rf <- Sys.time()

              ### trees ------
              use_trees <- if(isTRUE(trees)) {

                c(500, 1000, 2000, 4000)

              } else if(isFALSE(trees)) {

                500 # randomForest default for classification

              } else if(is.numeric(trees)) {

                trees

              }

              ### mtry ------
              use_mtry <- if(isTRUE(mtry)) {

                if(!prep$spatial_folds_used) {

                  # rf default for classification
                  seq(1, floor(sqrt(length(preds))), 1)

                } else {

                  # limit if spatial folds used
                  xm <- seq(1, floor(sqrt(length(preds))), 1)

                  xm <- xm[xm <= limit_spat_mtry]

                }

              } else if(isFALSE(mtry)) {

                floor(sqrt(length(preds))) # randomForest default for classification

              } else if(is.numeric(mtry)) {

                mtry

              }

              ### nodesize -----
              use_nodesize <- if(isTRUE(nodesize)) {

                seq(1, min(10, ceiling(min(start_df$n_p_train) / 2)), 2)

              } else if(isFALSE(nodesize)) {

                1 # randomForest default for classification

              } else if(is.numeric(nodesize)) {

                nodesize

              }

              safe_rf <- purrr::safely(randomForest::randomForest)

              tune_rf <- start_df %>%
                dplyr::cross_join(tibble::tibble(trees = use_trees)) %>%
                dplyr::cross_join(tibble::tibble(mtry = use_mtry)) %>%
                dplyr::cross_join(tibble::tibble(nodesize = use_nodesize)) %>%
                dplyr::mutate(tune_args = paste0("tr: ", trees
                                                 , ". mt: ", mtry
                                                 , ". ns: ", nodesize
                                                 )
                              , pa_train = purrr::map(pa_train, \(x) as.factor(x))
                              ) %>%
                dplyr::mutate(m = purrr::pmap(list(pa_train
                                                   , data_train
                                                   , n_p_train
                                                   , trees
                                                   , nodesize
                                                   , mtry
                                                   )
                                              , \(a, b, c, d, e, f) safe_rf(x = b
                                                                            , y = a
                                                                            , strata = a
                                                                            , sampsize = c(c, c)
                                                                            , ntree = d
                                                                            , nodesize = e
                                                                            , mtry = f
                                                                            )
                                              )
                              )

              probs <- tune_rf %>%
                dplyr::mutate(e = purrr::map(m, "error")) %>%
                dplyr::filter(purrr::map_lgl(e, \(x) !is.null(x)))

              if(nrow(probs)) {

                if(nrow(probs)) {

                purrr::map2(probs$e
                            , probs$tune_args
                            , \(a, b) {
                              readr::write_lines(paste0("error with tune_args "
                                                        , b
                                                        , ": "
                                                        , as.character(a)
                                                        )
                                                 , file = log_file
                                                 , append = TRUE
                                                 )
                              }
                           )

              }

              }

              tune$tune_rf <- tune_rf %>%
                dplyr::mutate(m = purrr::map(m, "result")) %>%
                dplyr::filter(purrr::map_lgl(m, \(x) !is.null(x))) %>%
                dplyr::mutate(e = purrr::pmap(list(m
                                                     , p_data_test
                                                     , a_data_test
                                                     )
                                                , \(a, b, c) evaluate_sdm(a, b, c
                                                                          , ...
                                                                          , do_gc = do_gc
                                                                          )
                                                )
                              ) %>%
                {if(keep_model) (.) %>% dplyr::select(! dplyr::where(is.list), m, e) else (.) %>% dplyr::select(! dplyr::where(is.list), e)}

              readr::write_lines(paste0("rf tune finished in: "
                                        , round(difftime(Sys.time(), start_rf, units = "mins"), 2)
                                        , " minutes"
                                        )
                                 , file = log_file
                                 , append = TRUE
                                 )

            }

          }

          if(best_run) {

            prep$blocks$block <- old_prep_block

          }

      } else {

        stop("no 'blocks' element in supplied prep")

      }

    } else {

      readr::write_lines("tune abandoned"
                         , file = log_file
                         , append = TRUE
                         )

    }


    # find best ------

    run <- all(!prep$abandoned
               , prep$finished
               , if(exists("tune_mean", where = tune, inherits = FALSE)) force_new else TRUE
               )

    if(run) {

      if(is.null(metrics_df)) stop("Can't find best model without metrics_df") else {

        # tunes ------
        tunes <- tune[grepl(paste0(algo, collapse = "|"), names(tune))] %>%
          dplyr::bind_rows(.id = "algo") %>%
          dplyr::mutate(algo = gsub("tune_", "", algo))

        if(nrow(tunes) > 0) {

          keeps <- c("algo", "spatial", "tune_args", "tunes"
                     , "fc", "rm", "treshold", "trees", "nodesize"
                     )

          metrics_df <- metrics_df %>%
            dplyr::mutate(summary_mets = metric %in% use_metrics)

          ## model stats-----
          stats <- if(any(!metrics_df$is_thresh[metrics_df$summary_mets])) {

            tunes %>%
              dplyr::mutate(stats = purrr::map(e
                                              , "stats"
                                              )
                            ) %>%
              tidyr::unnest(cols = c(stats)) %>%
              dplyr::group_by(dplyr::across(tidyselect::any_of(keeps))) %>%
              dplyr::mutate(k = as.factor(k)) %>%
              dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean)
                               , tunes = dplyr::n()
                               ) %>%
              dplyr::ungroup()

          } else tibble::tibble(algo = unique(tunes$algo))

          ## threshold stats-----
          tr_stats <- if(any(metrics_df$is_thresh[metrics_df$summary_mets])) {

            tunes %>%
              dplyr::mutate(stats = purrr::map(e
                                              , "tr_stats"
                                              )
                            ) %>%
              tidyr::unnest(cols = c(stats)) %>%
              dplyr::group_by(dplyr::across(tidyselect::any_of(keeps))) %>%
              dplyr::mutate(k = as.factor(k)) %>%
              dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean)) %>%
              dplyr::ungroup()

          } else tibble::tibble(algo = unique(tunes$algo))


          res_df <- tr_stats %>%
            dplyr::left_join(stats)


          if(nrow(res_df)) {

            ## thresholds----
            thresholds <- tunes %>%
              dplyr::mutate(stats = purrr::map(e
                                               , "thresholds"
                                               )
                            ) %>%
              tidyr::unnest(cols = c(stats)) %>%
              dplyr::group_by(dplyr::across(tidyselect::any_of(keeps))) %>%
              dplyr::mutate(k = as.factor(k)) %>%
              dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean)) %>%
              dplyr::ungroup()

            tune$tune_mean <- res_df %>%
              envFunc::make_metric_df(mets_df = metrics_df
                                      , context = keeps
                                      , mets_col = "summary_mets"
                                      , best_thresh = 1
                                      )  %>%
              dplyr::select(tidyselect::any_of(keeps), metric, value, combo, best) %>%
              dplyr::distinct() %>%
              tidyr::pivot_wider(names_from = "metric"
                                 , values_from = "value"
                                 ) %>%
              dplyr::arrange(desc(combo)) %>%
              dplyr::left_join(thresholds)

          }

        }

      }

    }

    # save -------
    # export before gc()
    tune$finished <- TRUE
    tune$log <- if(file.exists(log_file)) readr::read_lines(log_file) else NULL

    rio::export(tune, tune_file)

    # clean up --------

    if(do_gc) {

      stuff <- ls()

      delete_stuff <- stuff[! stuff %in% c(return_val, "return_val")]

      rm(list = delete_stuff)

      gc()

    }

    res <- if(return_val == "tune") get("tune") else list(tune_file = get("tune_file"))

    return(res)

  }




