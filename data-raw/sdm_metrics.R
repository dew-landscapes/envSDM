# focus on getting the presences right
  sdm_metrics <- tibble::tribble(
    ~metric, ~high_good, ~is_thresh
    # predicts::pa_evaluate
    , "TPR", T, T
    , "TNR", T, T
    , "FPR", F, T
    , "FNR", F, T
    , "PPP", T, T
    , "NPP", T, T
    , "MCR", F, T
    , "max_spec_sens", T, T
    , "no_omission", T, T
    , "equal_prevalence", T, T
    , "equal_sens_spec", T, T
    , "auc_po", T, F
    , "ODP", T, F
    # generated
    , "or10", T, T
    # flexsdm::fine_eval
    , "CBI", T, F
    , "CBI_rescale", T, F
    , "IMAE", T, F
    , "auc_po_flexsdm", T, F
    ) |>
    dplyr::mutate(within_mets = T
                  , across_mets = T
                  )
