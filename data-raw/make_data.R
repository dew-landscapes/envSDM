
  rm(list = ls())

  codes <- grep("make_data"
                , list.files(path = "data-raw",pattern = "\\.R$", full.names = TRUE)
                , value = TRUE
                , invert = TRUE
                )

  lapply(codes, source)

  datas <- ls(pattern = "sdm")

  do.call(save, c(lapply(datas,as.name), file = "data/data.rda"))


