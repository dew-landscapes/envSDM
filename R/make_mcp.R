
#' Minimum convex polygon around records
#'
#' For use as predict boundary for a taxa, but also see
#' `make_predict_boundary()`
#'
#' @param presence Cleaned and filtered dataframe of presences
#' @param out_file Character. Path for the mcp to be saved. Will be saved by
#' `sfarrow::st_write_parquet()`. Currently will not work very well with any
#' full stop in the path. Other file types are changed to .parquet
#' @param force_new Logical. If `out_file` exists, recreate it?
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates
#' @param in_crs epsg code for coordinates in `presence`
#' @param out_crs epsg code for coordinates in output mcp. Usually the same as
#' predictors
#' @param buf Distance to buffer the mcp. Passed to the `dist` argument of
#' `sf:st_buffer()` so units should follow that guidance.
#' @param clip sf to clip the mcp
#'
#' @return sf. `out_file` saved.
#' @export
#'
#' @example inst/examples/make_mcp_ex.R
  make_mcp <- function(presence
                       , out_file
                       , force_new = FALSE
                       , pres_x = "long"
                       , pres_y = "lat"
                       , in_crs = 4326
                       , out_crs = in_crs
                       , buf = 0
                       , clip = NULL
                       ) {

    run <- if(file.exists(out_file)) force_new else TRUE

    out_file <- gsub(tools::file_ext(out_file), "", out_file)
    out_file <- gsub("\\.$", "", out_file)
    out_file <- paste0(out_file, ".parquet")

    if(run) {

      fs::dir_create(dirname(out_file))

      res <- presence %>%
        dplyr::distinct(!!rlang::ensym(pres_y), !!rlang::ensym(pres_x)) %>%
        sf::st_as_sf(coords = c(pres_x, pres_y)
                     , crs = in_crs
                     ) %>%
        sf::st_transform(crs = out_crs) %>%
        sf::st_union() %>%
        sf::st_convex_hull() %>%
        sf::st_sf() %>%
        sf::st_buffer(buf)

      if(!is.null(clip)) {

        res <- res %>%
          sf::st_intersection(clip %>%
                                sf::st_transform(crs = out_crs) %>%
                                sf::st_make_valid()
                              )

      }

      if(isTRUE(nrow(res)>0)) sfarrow::st_write_parquet(res, out_file) else res <- NA

    } else {

      res <- sfarrow::st_read_parquet(out_file)

    }

    return(res)

  }
