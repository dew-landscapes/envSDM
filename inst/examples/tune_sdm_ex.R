
  # setup -------
  source(fs::path("inst", "examples", "prep_sdm_ex.R")) # need 'data' object

  purrr::map(data$out_dir
              , \(x) tune_sdm(prep = fs::path(x, "prep.rds")
                              , out_dir = x
                              , fc = "lq"
                              , rm = c(1, 2)
                              , trees = 500
                              , mtry = c(1:3)
                              , nodesize = c(1, 3)
                              , limit_p = 3
                              , use_metrics = c("auc_po", "CBI_rescale", "IMAE")
                              #, force_new = TRUE
                              )
              )

  # which tune args were best using 'combo'?
  # BUT, not a sensible comparison as, between rows, the models are not all built on the same data!
  data %>%
    dplyr::mutate(tune = purrr::map(tune, rio::import, trust = TRUE)
                  , tune_mean = purrr::map(tune, "tune_mean")
                  ) %>%
    tidyr::unnest(cols = c(tune_mean)) %>%
    dplyr::filter(best) |> # used 'combo' to determine 'best' as default in tune_sdm
    dplyr::select(taxa, algo, tidyselect::any_of(names(sdms)), tidyselect::where(is.numeric))

  # random forest nearly always 'won' with both nodeside and mtry varying
