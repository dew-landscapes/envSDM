
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  data <- fs::path(system.file(package = "envSDM"), "examples") |>
    fs::dir_ls(regexp = "prep\\.rds$"
               , recurse = TRUE
               ) |>
    tibble::enframe(name = NULL, value = "prep") |>
    dplyr::mutate(taxa = gsub("\\.rds", "", basename(dirname(prep)))
                  , out_dir = fs::path(out_dir, taxa)
                  )

  purrr::map(data$out_dir
              , \(x) tune_sdm(prep = fs::path(x, "prep.rds")
                              , out_dir = x
                              , fc = "lq"
                              , rm = c(2, 3)
                              , trees = 500
                              , mtry = c(1:3)
                              , nodesize = c(1, 2, 3)
                              , limit_p = 3
                              , use_metrics = c("auc_po", "CBI_rescale", "IMAE", "or10")
                              #, force_new = TRUE
                              )
              )

  # which tune args were best for each taxa using 'combo'?
  data %>%
    dplyr::mutate(tune = fs::path(out_dir, "tune.rds")
                  , tune = purrr::map(tune, rio::import, trust = TRUE)
                  , tune_mean = purrr::map(tune, "tune_mean")
                  ) %>%
    tidyr::unnest(cols = c(tune_mean)) %>%
    dplyr::filter(best) %>% # used 'combo' to determine 'best' as default in tune_sdm
    dplyr::select(taxa, algo, tune_args, combo, auc_po, IMAE, CBI, max_spec_sens)

  # or best tune args choosing on just auc_po?
  data %>%
    dplyr::mutate(tune = fs::path(out_dir, "tune.rds")
                  , tune = purrr::map(tune, rio::import, trust = TRUE)
                  , all = purrr::map(tune, "tune_mean")
                  ) %>%
    tidyr::unnest(cols = c(all)) %>%
    dplyr::group_by(taxa) %>%
    dplyr::filter(auc_po == max(auc_po)) %>%
    dplyr::ungroup() %>%
    dplyr::select(taxa, algo, tune_args, auc_po, IMAE, CBI, max_spec_sens)
