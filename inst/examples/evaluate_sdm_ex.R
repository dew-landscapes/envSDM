
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  preps <- fs::dir_ls(out_dir, regexp = "prep.rds", recurse = TRUE)

  prep <- rio::import(preps[[1]], trust = TRUE)

  full_run <- rio::import(fs::path(dirname(preps[[1]]), "combo", "full_run.rds"), trust = TRUE)
  algo <- full_run$tune_mean$algo[[1]]
  model <- full_run[[paste0("tune_", algo)]]$m[[1]]

  presences <- prep$testing[prep$testing$pa == 1, ]
  background <- prep$testing[prep$testing$pa == 0, ]

  evaluate_sdm(full_run$tune_rf$m[[1]]
               , p_test = presences
               , b_test = background
               )
