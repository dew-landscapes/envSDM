
#' Merge polygons (or polygon files) to form a single minimum convex polygon
#' (mcp)
#'
#' Primary use is to create a predict boundary for `prep_sdm()`, merging (an
#' existing) minimum convex polygon around points with other sources of taxa
#' distribution. Optionally, applying a buffer around the resulting mcp; and
#' clipping to a (usually coastal) boundary. The predict boundary is then used
#' for generation of background points and for masking the 'full' predict (to
#' the full extent of the environmental variables.
#'
#'
#' @param poly_list List of paths or list of sf
#' @param out_file Character name of file to save
#' @param buffer_metres Numeric. Distance in metres to buffer the mcp
#' @param col_name Name of column to create in the resulting mcp
#' @param col_name_val Value to provide in the column in the resulting mcp
#' @param clip sf. Clip the resulting mcp back to this.
#' @param out_crs Numeric. [epsg](https://epsg.io/) code
#' @param return_poly Logical. Return the mcp, or alternatively `invisible(NULL)`
#' @param force_new Logical. If `out_file` exists, recreate it?
#'
#' @return `if(return_poly) sf else invisible(NULL)`. .parquet mcp written to
#' `out_file`
#' @export
#'
#' @examples
make_predict_boundary <- function(poly_list
                                  , out_file
                                  , buffer_metres = 0
                                  , col_name = "taxa"
                                  , col_name_val = basename(dirname(out_file))
                                  , clip = NULL
                                  , out_crs
                                  , return_poly = FALSE
                                  , force_new = FALSE
                                  ) {

  run <- if(base::file.exists(out_file)) force_new else TRUE

  if(run) {

    if(all(purrr::map_lgl(poly_list, is.character))) {

      poly_list <- stats::na.omit(unlist(poly_list))
      poly_list <- poly_list[base::file.exists(poly_list)]

      polys <- purrr::map(poly_list
                          , \(x) sfarrow::st_read_parquet(x)
                          )

      polys <- polys[purrr::map_lgl(polys, \(x) nrow(x) > 0)]

    } else {

      polys <- poly_list

    }

    if(length(polys)) {

      fs::dir_create(dirname(out_file))

      mcp <- polys %>%
        purrr::map(\(x) sf::st_geometry(x) %>%
                     sf::st_as_sf() %>%
                     sf::st_transform(crs = out_crs)
                   ) %>%
        dplyr::bind_rows() %>%
        sf::st_cast("MULTIPOINT") %>%
        dplyr::summarise() %>%
        sf::st_convex_hull() %>%
        sf::st_buffer(buffer_metres) %>%
        dplyr::mutate(!!rlang::ensym(col_name) := col_name_val)

      if(!is.null(clip)) {

        clip <- clip %>%
          sf::st_transform(crs = out_crs)

        mcp <- sf::st_intersection(mcp, clip)


      }

      sfarrow::st_write_parquet(mcp, out_file)

    } else poly <- NA

  } else {

    if(return_poly) mcp <- sfarrow::st_read_parquet(out_file)

  }

  if(return_poly) return(mcp) else return(invisible(NULL))

}
