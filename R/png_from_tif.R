
  png_from_tif <- function(x # raster or raster path
                           , title = NULL
                           , dots = NULL
                           , trim = TRUE
                           , out_png = NULL
                           , do_gc = FALSE
                           ) {

    if(is.character(x)) {

      x <- terra::rast(x)

    }

    if(is.null(out_png)) out_png <- gsub("\\..*$", ".png", terra::sources(x))

    png(filename = out_png
        , width = 8
        , height = 8
        , units = "in"
        , bg = "white"
        , res = 300
        )

    safe_trim <- purrr::safely(terra::trim)

    r <- if(trim) safe_trim(x) else list(result = x)

    r <- if(is.null(r$error)) {

      r$result

    } else x

    terra::plot(r
                , main = title
                )

    if("sf" %in% class(dots)) {

      terra::plot(terra::vect(dots)
                  , pch = 4
                  , cex = 0.5
                  , col = scales::alpha("blue", 0.5)
                  , add = TRUE
                  )

    }

    dev.off()

    if(do_gc) {

      rm(list = ls())

      gc()

    }

  }
