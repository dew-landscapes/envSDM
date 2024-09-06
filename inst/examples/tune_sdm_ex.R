
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  data <- file.path(system.file(package = "predicts"), "ex") |>
    fs::dir_ls(regexp = "\\.csv$") |>
    tibble::enframe(name = NULL, value = "path") |>
    dplyr::mutate(taxa = gsub("\\.csv", "", basename(path))
                  , presence = purrr::map(path, rio::import, setclass = "tibble")
                  , presence = purrr::map(presence
                                          , \(x) x |>
                                            dplyr::filter(!is.na(lat)
                                                          , !is.na(lon)
                                                          )
                                          )
                  , out_dir = fs::path(out_dir, taxa)
                  )

  purrr::walk(data$out_dir
              , \(x) tune_sdm(out_dir = x)
              )

  # which tune args were best for each taxa using 'combo'?
  data %>%
    dplyr::mutate(eval_file = fs::path(out_dir, "evaluation.csv")
                  , eval = purrr::map(eval_file, rio::import, setclass = "tibble")
                  ) %>%
    tidyr::unnest(cols = c(eval)) %>%
    dplyr::filter(best) %>%
    dplyr::select(taxa, algo, tune_args, combo, auc_po, IMAE, CBI, max_spec_sens)

  # or best tune args choosing on just auc_po?
  data %>%
    dplyr::mutate(eval_file = fs::path(out_dir, "evaluation.csv")
                  , eval = purrr::map(eval_file, rio::import, setclass = "tibble")
                  ) %>%
    tidyr::unnest(cols = c(eval)) %>%
    dplyr::group_by(taxa) %>%
    dplyr::filter(auc_po == max(auc_po)) %>%
    dplyr::ungroup() %>%
    dplyr::select(taxa, algo, tune_args, combo, auc_po, IMAE, CBI, max_spec_sens)
