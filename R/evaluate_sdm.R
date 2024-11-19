

#' Evaluate an SDM
#'
#' Returns various evaluation metrics from `predicts::pa_evaluate()` and
#' `flexsdm::sdm_eval()`.
#'
#' @param m SDM result within `tune_sdm()`
#' @param p_test Presence test data generated within `tune_sdm()`
#' @param b_test Background test data generated within `tune_sdm()`
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel. Note, actually usees `rm(list = ls(pattern = "^[^e$]"))`.
#' @param ... Passed to both `terra::predict()` and `predicts::pa_evaluate()`
#'
#' @return paModelEvaluation (see `predicts::pa_evaluate()`) with extra metrics
#' from `flexsdm::sdm_eval()`: AUC as auc_po_flexsdm; BOYCE as CBI; CBI_rescale
#' (CBI is -1 to 1, CBI_rescale is 0 to 1); and IMAE.
#' @export
#' @keywords internal
#'
#' @example inst/examples/evaluate_sdm_ex.R
  evaluate_sdm <- function(m
                           , p_test # presence test data
                           , b_test # background (or absence) test data
                           , do_gc = FALSE
                           , ...
                           ) {

    if(any(class(m) %in% c("maxnet", "envelope_model"))) {

      p_type <- "cloglog"
      index <- 1

    } else {

      p_type <- "prob"
      index <- "1"

    }

    # predict not working for envelope models. Don't think it could find the appropriate predict. no issue for rf or maxnet
    requireNamespace("predicts", quietly = TRUE)

    # requireNamespace not fixing issue
    if("envelope_model" %in% class(m)) library("predicts")

    new_data <- dplyr::bind_rows(p_test, b_test)

    p <- predict(m
                 , type = p_type
                 , newdata = new_data
                 , x = new_data
                 , ...
                 ) %>%
      tibble::as_tibble(.name_repair = "minimal") %>%
      dplyr::select(pred = tidyselect::any_of(index)) %>%
      dplyr::mutate(p = c(rep(1, nrow(p_test)), rep(0, nrow(b_test))))

    e <- predicts::pa_evaluate(p[p$p == 1, 1][[1]]
                               , p[p$p == 0, 1][[1]]
                               , type = p_type
                               , ...
                               )

    e_fsdm <- flexsdm::sdm_eval(p[p$p == 1, 1][[1]]
                               , p[p$p == 0, 1][[1]]
                               )

    e@stats$auc_po <- e@stats$auc
    e@stats$auc_po_flexsdm <- mean(e_fsdm$AUC, na.rm = TRUE)
    e@stats$CBI <- mean(e_fsdm$BOYCE, na.rm = TRUE)
    e@stats$CBI_rescale <- (e@stats$CBI + 1) / 2
    e@stats$IMAE <- mean(e_fsdm$IMAE, na.rm = TRUE)

    if(do_gc) {

      rm(list = ls(pattern = "^[^e$]"))

      gc()

    }

    return(e)

  }

