
  if(FALSE) {

    # clean up (delete) examples
    fs::dir_info("inst/examples") |>
      dplyr::filter(type == "directory") |>
      dplyr::pull(path) |>
      grep("tif$", x = _, value = TRUE, invert = TRUE) |>
      fs::dir_delete()

  }

  devtools::load_all()
  source("inst/examples/prep_sdm_ex.R")
  source("inst/examples/tune_sdm_ex.R")
  source("inst/examples/predict_sdm_ex.R")
  source("inst/examples/thresh_sdm_ex.R")
  source("inst/examples/evaluate_sdm_ex.R")
