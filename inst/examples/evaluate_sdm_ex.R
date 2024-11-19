
  out_dir <- file.path(system.file(package = "envSDM"), "examples", "acaule")

  prep <- rio::import(fs::path(out_dir, "prep.rds"))

  presences <- prep$blocks[prep$blocks$pa == 1 & prep$blocks$block != 1, ]
  background <- prep$blocks[prep$blocks$pa == 0 & prep$blocks$block != 1, ]

  evaluate_ex_dir <- fs::path(out_dir, "evaluate_example")
  fs::dir_create(evaluate_ex_dir)

  tune_sdm(out_dir = out_dir
           , algo = "rf"
           , trees = 500
           , mtry = FALSE
           , nodesize = FALSE
           , keep_model = TRUE
           , save_to = evaluate_ex_dir
           )

  tune <- rio::import(fs::path(evaluate_ex_dir, "tune.rds")
                      , setclass = "tibble"
                      )

  evaluate_sdm(tune$m[[1]]
               , p_test = presences
               , b_test = background
               )
