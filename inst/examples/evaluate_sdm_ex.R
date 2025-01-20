
  out_dir <- file.path(system.file(package = "envSDM"), "examples")

  source(fs::path(out_dir, "tune_sdm_ex.R")) # make sure following prep file exists

  prep <- rio::import(fs::path(out_dir, "acaule", "prep.rds"))

  model <- tune_sdm(prep = prep
                    , out_dir = FALSE
                    , return_val = "object"
                    , algo = "rf"
                    , trees = 500
                    , mtry = 2
                    , nodesize = 1
                    , keep_model = TRUE
                    )

  presences <- prep$blocks[prep$blocks$pa == 1 & prep$blocks$block != 1, ]
  background <- prep$blocks[prep$blocks$pa == 0 & prep$blocks$block != 1, ]

  evaluate_sdm(model$tune_rf$m[[1]]
               , p_test = presences
               , b_test = background
               )
