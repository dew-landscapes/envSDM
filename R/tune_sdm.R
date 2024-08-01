
#' Tune, and evaluate, species distribution models
#'
#' @param this_taxa Character. Name of taxa. Used to name outputs. If `NULL`,
#' this will be `basename(dirname(out_dir))`.
#' @param out_dir Character. Name of directory into which results will be saved.
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
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param force_new Logical. If outputs already exist, should they be remade?
#' @param save_to Character. Name of path to save results. Defaults to out_dir
#' which works if not `best_run`. For a best run from `run_full_sdm()`, this
#' will default to `fs::path(out_dir, "best")`. Otherwise set as anything.
#' @param ... Passed to `evaluate_sdm()`
#'
#' @return
#' @export
#'
#' @examples man/examples/tune_sdm_ex.R
  tune_sdm <- function(this_taxa = NULL
                       , out_dir
                       , algo = c("all", "maxnet", "bioclim", "envelope", "rf") # envelope = bioclim
                       , fc = "auto_feature" # maxnet tune
                       , limit_p = FALSE # T, F or number of preds above which to limit p
                       , rm = seq(1, 6, 0.5) # maxnet tune
                       , trees = c(500, 1000, 2000) # T, F or numeric
                       , mtry = TRUE # T, F or numeric
                       , nodesize = c(1, 2) # T, F or numeric
                       , keep_model = FALSE
                       , best_run = FALSE
                       , metrics_df = envSDM::sdm_metrics
                       , do_gc = TRUE
                       , force_new = FALSE
                       , save_to = out_dir
                       , ...
                       ) {

    # files-----
    ## already done
    prep_log <- fs::path(out_dir
                           , "prep.log"
                           )

    prep_file <- fs::path(out_dir, "prep.rds")

    ## new
    tune_file <- fs::path(save_to, "tune.rds")

    eval_file <- fs::path(out_dir
                          , "evaluation.csv"
                          ) # eval_file is always in taxa folder, not save_to folder

    tune_log <- fs::path(save_to
                         , "tune.log"
                         )

    # run?-----
    prep_log_text <- paste0(readLines(prep_log), collapse = " ")

    run <- all(file.exists(prep_file) # need prep or can't run
               , !grepl("SDM abandoned", prep_log_text) # too few records can't run
               , grepl("prep end", prep_log_text) # prep needs to have finished
               )

    if(run) {

      if(!best_run) message("tuning ", basename(out_dir))

      run_tune <- if(file.exists(tune_file)) force_new else TRUE

      run_find_best <- if(file.exists(eval_file)) force_new else TRUE

      if(any(run_tune, run_find_best)) {

        prep <- rio::import(prep_file)

      }

      # tune -----
      if(run_tune) {

        if(exists("blocks", where = prep)) {

          tunes <- list()

          # start timer ------
          tune_timer <- envFunc::timer(process = "tune start"
                              , file = if(best_run) "best" else "tune"
                              , time_df = NULL
                              , log = tune_log
                              , write_log = TRUE
                              )


          ## setup ------

          if(is.null(this_taxa)) this_taxa <- basename(dirname(prep_file))

          nobs <- nrow(prep$blocks[prep$blocks$pa == 1,])

          ## start data frame -----

          if(best_run) {

            old_prep_block <- prep$blocks$block

            prep$blocks$block <- 1

          }

          single_block <- if(length(unique(prep$blocks[prep$blocks$pa == 1,]$block)) == 1) TRUE else FALSE

          preds <- prep$correlated$env_cols[!prep$correlated$env_cols %in% prep$correlated$remove_env]

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
                          , n_p_train = purrr::map_dbl(pa_train, sum)
                          , n_a_train = purrr::map2_dbl(data_train, n_p_train
                                                        , \(x, y) nrow(x) - y
                                                        )
                          # 'hack' to ensure min of p's and a's is used (as n_p_train)
                          , n_p_train = purrr::map2_dbl(n_p_train
                                                        , n_a_train
                                                        , \(x, y) min(x, y)
                                                        )
                          ) %>%
            dplyr::filter(n_p_test > 0)

          if(nrow(start_df) < length(unique(prep$blocks$block))) {

            # This should now be handled within prep, so should not be needed here

            tune_timer <- envFunc::timer("warning"
                                , notes = paste0("Some folds had zero presences. These were removed, leaving only ", nrow(start_df), " folds")
                                , time_df = tune_timer
                                )

          }

          ## tune maxnet------

          if(any(c("all", "maxnet") %in% algo)) {

            message("maxnet tune")

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

            # if limit_p, ensure p not in fc
            if(!isFALSE(limit_p)) {

              # Only limit if too many preds, or limit_p is TRUE
              if(any(length(preds) > limit_p, isTRUE(limit_p))) {
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
                            tune_timer <- envFunc::timer("error"
                                                , notes = paste0("tune_args: ", b, ": ", a)
                                                , time_df = tune_timer
                                                , write_log = TRUE
                                                )
                            }
                         )

            }

            tunes$tune_maxnet <- tune_maxnet %>%
              dplyr::mutate(m = purrr::map(m, "result")) %>%
              dplyr::filter(purrr::map_lgl(m, \(x) !is.null(x))) %>%
              dplyr::mutate(e = purrr::pmap(list(m
                                                   , p_data_test
                                                   , a_data_test
                                                   )
                                              , evaluate_sdm
                                              , ...
                                              , do_gc = do_gc
                                              )
                            ) %>%
              {if(keep_model) (.) %>% dplyr::select(! dplyr::where(is.list), m, e) else (.) %>% dplyr::select(! dplyr::where(is.list), e)}

            tune_timer <- envFunc::timer("maxnet"
                                , time_df = tune_timer
                                , write_log = TRUE
                                )

          }


          ## tune envelope -------

          if(any(c("all", "envelope", "bioclim", "env") %in% algo)) {

            message("envelope tune")

            tunes$tune_env <- start_df %>%
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

            tune_timer <- envFunc::timer("envelope"
                                , time_df = tune_timer
                                , write_log = TRUE
                                )

          }

          ## tune rf --------

          if(any(c("all", "rf", "randomForest") %in% algo)) {

            message("rf tune")

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

                xm <- xm[xm <= 3]

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

              purrr::map2(probs$e
                          , probs$tune_args
                          , \(a, b) {
                            tune_timer <- envFunc::timer("error"
                                                , notes = paste0("tune_args: ", b, ": ", a)
                                                , time_df = tune_timer
                                                , write_log = TRUE
                                                )
                            }
                         )

            }

            tunes$tune_rf <- tune_rf %>%
              dplyr::mutate(m = purrr::map(m, "result")) %>%
              dplyr::filter(purrr::map_lgl(m, \(x) !is.null(x))) %>%
              dplyr::mutate(e = purrr::pmap(list(m
                                                 , p_data_test
                                                 , a_data_test
                                                 )
                                            , evaluate_sdm
                                            , ...
                                            , do_gc = do_gc
                                            )
                            ) %>%
              {if(keep_model) (.) %>% dplyr::select(! dplyr::where(is.list), m, e) else (.) %>% dplyr::select(! dplyr::where(is.list), e)}

            tune_timer <- envFunc::timer("rf"
                                , time_df = tune_timer
                                , write_log = TRUE
                                )

          }

          if(best_run) {

            prep$blocks$block <- old_prep_block

          }

          ## tunes ------

          tunes <- tunes %>%
            dplyr::bind_rows(.id = "algo") %>%
            dplyr::mutate(algo = gsub("tune_", "", algo))

          rio::export(tunes
                      , tune_file
                      )

        }

        # end timer ------
        if(exists("tune_timer")) {

          tune_timer <- envFunc::timer("tune end"
                                , time_df = tune_timer
                                )

        }

      }

      if(run_find_best) {

        # find best ------

        message(paste0("find best model for ", basename(out_dir)))

        run <- if(file.exists(eval_file)) force_new else TRUE

        if(run) {

          if(!exists("tunes", where = environment(), inherits = FALSE)) {

            tunes <- rio::import(fs::path(out_dir, "tune.rds"))

          }

          if(is.null(metrics_df)) warning("Can't find best model without metrics") else {

            if(nrow(tunes) > 0) {

              keeps <- c("algo", "spatial", "tune_args", "tunes"
                         , "fc", "rm", "treshold", "trees", "nodesize"
                         )

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

                best <- res_df %>%
                  envFunc::make_metric_df(mets_df = metrics_df
                                          , context = keeps
                                          , mets_col = "summary_mets"
                                          , best_thresh = 1
                                          ) %>%
                  dplyr::select(tidyselect::any_of(keeps), metric, value, combo, best) %>%
                  dplyr::distinct() %>%
                  tidyr::pivot_wider(names_from = "metric"
                                     , values_from = "value"
                                     ) %>%
                  dplyr::arrange(desc(combo)) %>%
                  dplyr::left_join(thresholds)

                best_tune <- res_df %>%
                  dplyr::inner_join(best %>%
                                      dplyr::select(tidyselect::any_of(keeps)
                                                    , tidyselect::any_of(names(thresholds))
                                                    , best
                                                    ) %>%
                                      dplyr::filter(best)
                                    )

                rio::export(res_df %>%
                              dplyr::select(! dplyr::where(is.list))
                            , fs::path(save_to
                                       , "metrics.csv"
                                       )
                            )

                rio::export(best
                            , eval_file
                            )

                tune_timer <- envFunc::timer("find best"
                                    , time_df = tune_timer
                                    , write_log = TRUE
                                    )

              }

            }

          }

        }

      }

      # clean up --------

      if(do_gc) {

        rm(list = ls())

        gc()

      }

    }

  return(invisible(NULL))

  }




