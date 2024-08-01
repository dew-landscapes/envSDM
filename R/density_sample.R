
#' Create a spatial sample
#'
#' Same idea as `terra::spatSample(weights = weight_raster)`, but can use a
#' density raster with a different spatial resolution. Can be much quicker than
#' `terra::spatSample()` for large areas and small resolutions.
#'
#' @param x spatRaster. Sample will be cell XY values from this raster.
#' @param dens_rast spatRaster. Density raster. This is equivalent to the
#' `weights` argument of `terra::spatSample()` but does not have to be at the
#' same resolution as `x`. Usually made within `prep_sdm()`
#' @param samp_df Dataframe. Usually built by `prep_sdm()`
#' @param max_sample Numeric. How many iterations to try to meet the density
#' targets in `samp_df`.
#' @param verbose Logical. Prints how many iterations have been attempted.
#'
#' @return Dataframe
#' @export
#' @keywords internal
#'
#' @examples
  density_sample <- function(x
                             , dens_rast
                             , samp_df
                             , max_sample = 200
                             , verbose = FALSE
                             ) {

    a_sample <- function(x, dens_rast, samp_df) {

      sample(1:(terra::ncell(x)), sum(samp_df$target)) %>%
        terra::xyFromCell(x, .) %>%
        as.data.frame %>%
        cbind(terra::extract(dens_rast, as.matrix(.))) %>%
        stats::na.omit() %>%
        dplyr::rename(value = 3) %>%
        merge(samp_df[, c("value", "target")])

    }

    new <- a_sample(x, dens_rast, samp_df)

    todo <- TRUE
    counter <- 1

    while(all(counter < max_sample, todo)) {

      new <- rbind(new[,1:4], a_sample(x, dens_rast, samp_df)) %>%
        dplyr::add_count(value
                         , name = "n"
                         ) %>%
        dplyr::group_by(value) %>%
        dplyr::slice(1:unique(target)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct() %>%
        dplyr::add_count(value, name = "n")

      todo <- new %>%
        dplyr::count(value, target) %>%
        dplyr::mutate(todo = n < target) %>%
        dplyr::pull(todo) %>%
        sum()

      counter <- counter + 1

      if(verbose) message("counter: ", counter)

    }

    return(new)

  }
