
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
