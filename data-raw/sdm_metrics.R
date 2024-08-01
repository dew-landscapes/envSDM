# focus on getting the presences right
  sdm_metrics <- tibble::tibble(metric = c(
  # predicts::pa_evaluate
    "TPR", "TNR"
    , "FPR", "FNR"
    , "PPP", "NPP"
    , "MCR", "auc_po"
    , "ODP", "kappa"
  # flexsdm::fine_eval
    , "CBI_rescale", "IMAE"
    , "auc_po_flexsdm", "CBI"
    )
    , high_good = c(T, T
                    , F, F
                    , T, T
                    , F, T
                    , T, T
                    , T, T
                    , T, T
                    )
  , is_thresh = c(T, T
                  , T, T
                  , T, T
                  , F, F
                  , F, F
                  , F, F
                  , F, F
                  )
    , summary_mets = c(F, F
                       , F, F
                       , F, F
                       , F, T
                       , F, F
                       , T, T
                       , F, F
                       )
    , within_mets = T
    , across_mets = T
    )
