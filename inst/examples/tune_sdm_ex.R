
out_dir <- file.path(system.file(package = "envSDM"), "examples")

# setup -------
data <- readRDS(fs::path(out_dir, "data.rds"))

max_cores <- nrow(data)
use_cores <- min(max_cores, parallel::detectCores() - 1)

future::plan(future::multisession(workers = use_cores))

furrr::future_walk(data$out_dir
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

future::plan(future::sequential())

# which tune args were 'best' using combo?
# BUT, possibly spurious comparison as, between rows, the models are not all built on the same data! Same presences though.
data |>
  dplyr::mutate(tuned = purrr::map_lgl(tune, file.exists)) |>
  dplyr::filter(tuned) |>
  dplyr::mutate(tune = purrr::map(tune, rio::import, trust = TRUE)
                , tune_mean = purrr::map(tune, "tune_mean")
                ) |>
  tidyr::unnest(cols = c(tune_mean)) |>
  dplyr::filter(combo == max(combo), .by = taxa) |> # used 'combo' to determine 'best' as default in tune_sdm
  dplyr::select(taxa, algo, hold_prop, stretch, spatial_folds, tidyselect::where(is.numeric))
