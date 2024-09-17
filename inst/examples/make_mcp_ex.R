
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
                , out_mcp = fs::path(out_dir, "mcp.parquet")
  )


# mcps --------

purrr::pwalk(list(data$presence
                  , data$out_mcp
)
, \(x, y) make_mcp(x, y, pres_x = "lon")
)
