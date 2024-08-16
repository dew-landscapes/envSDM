
#' Create a .png from a .tif
#'
#'
#'
#' @param x spatRaster or path to .tif
#' @param title Character. Title to add to the .png
#' @param dots sf. Usually presences. Added as points.
#' @param trim Logical. Run `terra::trim()` before writing to .png?
#' @param out_png Character. Name of .png file to save. If `NULL` will be the
#' same file name as `terra::sources(x)` with the file type as .png
#' @param do_gc  Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#'
#' @return
#' @export
#'
#' @examples
  png_from_tif <- function(x
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

    return(invisible(NULL))

  }
