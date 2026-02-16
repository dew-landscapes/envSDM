
#' Prepare for running an SDM
#'
#' The background sampling includes code based on a
#' [Geographic Information Systems stack exchange](https://gis.stackexchange.com/)
#' [answer](https://gis.stackexchange.com/a/224347)
#' by user [Spacedman](https://gis.stackexchange.com/users/865/spacedman).
#'
#' If memory is an issue, try adjusting `terra_options` and/or `do_gc`.
#'
#' To help build the density raster for assigning background points, 'absence'
#' data can be supplied in `presence` as `0` values. e.g. For a bird, absence
#' data might be generated from other sites where birds were recorded but
#' `this_taxa` was not. Note that these `0` 'records' are not used directly as
#' background points but instead are used to improve the sampling density raster
#' against which background points are assigned.
#'
#' @param this_taxa Character. Name of taxa. Only used to print some messages.
#' Ignored if NULL
#' @param out_dir FALSE or character. If FALSE the result of prep_sdm will be
#' saved to a temporary folder. If character, a file 'prep.rds' will be created
#' at the path defined by out_dir.
#' @param return_val Character: "object" or "path". Both return a named list. In
#' the case of "path" the named list is simply list(prep = out_dir). Will be set
#' to "object" if `out_dir` is FALSE.
#' @param presence Dataframe of presences with columns `pres_x` and `pres_y`.
#' @param pres_col Character. Name of column in `presence` that defines presence
#' (`1`) or absence (`0`). Optional if only presence data is supplied.
#' @param pres_val Numeric. Values in `pres_col` that represent presences.
#' Optional if only presence data is supplied.
#' @param pres_crs Anything that will return a legitimate crs when passed to the
#' crs attribute of `sf::st_transform()` or `sf::st_as_sf()`.
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates
#' @param pred_limit Limit the background points and predictions?
#' Can be `TRUE` (use `presence` to generate a minimum convex polygon to use as
#' a limit. Not recommended as the points in `presence` have usually been
#' filtered to very accurate spatial reliability and thus may be missing a large
#' number of legitimate records); `FALSE` (the full extent of the predictors
#'  will be used); path to existing .parquet to use; or sf object.
#' @param limit_buffer Numeric. Apply this buffer to `pred_limit`. Only used if
#' `pred_limit` is `TRUE`. Passed to the `dist` argument of `sf::st_buffer()`.
#' @param pred_clip sf. Optional sf to clip the pred_limit back to (e.g. to
#' prevent prediction into ocean).
#' @param predictors Character. Vector of paths to predictor `.tif` files.
#' @param is_env_pred Logical. Does the naming of the directory and files in
#' `predictors` follow the pattern required by `envRaster::parse_env_tif()`?
#' @param cat_preds Character. Vector of predictor names that are character.
#' @param terra_options Passed to `terra::terraOptions()`. e.g. list(memfrac =
#' 0.6)
#' @param num_bg Numeric. How many background points?
#' @param prop_abs Character. Is `num_bg` a proportion (`prop`) of the number of
#' records in `presence` or an absolute (`abs`) number?
#' @param many_p_prop Numeric. Ensure the number of background points is at
#' least `many_p_prop * number of presences`. e.g. If there are more than 5000
#' presences and num_bg is set at `10000` and `many_p_prop` is `2`, then num_bg
#' will be increased to `many_p_prop * nrow(presences)`
#' @param folds Numeric. How many folds to use in cross validation? Will be
#' adjusted downwards if number of presences do not support `folds * min_fold_n`
#' @param spatial_folds Logical. Use spatial folds? Even if `TRUE`, can resort
#' to non-spatial cv if presences per fold do not meet `min_fold_n` or there are
#' not enough presences to support more than one fold.
#' @param repeats Numeric. Number of repeated cross validations.
#' @param folds_div Numeric. The square root of the predict area is divided
#' by this value before being passed to the `block_dist` argument of
#' `blockCV::cv_spatial()`. If using cross validation, `fold_div`
#' must be of the same length as `1:folds`.
#' @param max_repeat_corr Numeric. Maximum correlation allowed between the folds
#' of any two spatial folds before one of the correlated folds will be set to
#' non-spatial and the folds reallocated. Correlation tested on presences only.
#' @param min_fold_n Numeric. Sets both minimum number of presences, and,
#' by default, the minimum number of presences required for a model.
#' @param stretch_value Numeric. Stretch the density raster to this value.
#' @param dens_res `NULL` or numeric. Resolution (in metres) of density raster.
#' Set to `NULL` to use the same resolution as the predictors.
#' @param reduce_env_thresh_corr Numeric. Threshold used to flag highly
#' correlated variables. Set to 1 to skip this step. If > 0, highly
#' correlated and low importance variables will be removed. In the case of
#' highly correlated pairs of variables, only one is removed.
#' @param reduce_env_quant_rf_imp Numeric. Bottom quantile of importance values
#' to drop.
#' @param hold_prop Numeric. Proportion of data to be held back from training
#' to use to validate the final model.
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param force_new Logical. If outputs already exist, should they be remade?
#'
#' @return If `return_val` is "object" a named list. If `return_val` is "path"
#' a path to the saved file. If `out_dir` is a valid path, the 'full
#' result' (irrespective of `return_val`) is also saved to
#' `fs::path(out_dir, "prep.rds")`. The 'full result' is a named list with
#' elements:
#' * log:
#'     + a log of (rough) timings and other information from the process
#' * abandoned:
#'     + Logical indicating if the sdm was abandoned. If abandoned is TRUE, some
#'     list elements may not be present
#' * presence_ras:
#'     + tibble with two columns ('x' and 'y') representing unique cell
#'     centroids on the predictors at presences supplied in argument `presence`
#' * predict_boundary:
#'     + sf used to limit the background points and used by `predict_sdm()` to
#'     generate the 'mask'ed output
#'     + i.e. the predict boundary is both callibration and predict boundary
#' * bg_points:
#'    + sf of cell centroids representing unique cell centroids for background
#'    points
#' * folds
#'     + data.frame with columns:
#'       + `pa`: presence (1) or absence/background (0)
#'       + `x` and `y`: cell centroids for each presence and absence
#'       + `fold`: the spatial fold to which the row belongs
#'       +  a column with values for each of `predictors` at `x` and `y`
#' * spatial_folds_used:
#'     + logical indicating if spatial folds were used. This may differ from
#'     the `spatial_folds` argument provided to `prep_sdm()` if an attempt to
#'     use spatial folds failed to meet desired `folds` and `min_fold_n`
#' * correlated:
#'     + list with elements as per `envModel::reduce_env()`, or, if
#'     `reduce_env` is `FALSE`, a list with elements `remove_env` which is
#'     empty, and `env_var` and `keep`, which both contain the names of all
#'     predictors.
#'
#' @export
#'
#' @example inst/examples/prep_sdm_ex.R
  prep_sdm <- function(this_taxa = NULL
                       , out_dir = FALSE
                       , return_val = "path"
                       , presence
                       , pres_col = "pa"
                       , pres_val = 1
                       , pres_crs = 4326
                       , pres_x = "long"
                       , pres_y = "lat"
                       , pred_limit = TRUE
                       , limit_buffer = 0
                       , pred_clip = NULL
                       , predictors
                       , is_env_pred = TRUE
                       , terra_options = NULL
                       , cat_preds = NULL
                       , num_bg = 10000
                       , prop_abs = "abs"
                       , many_p_prop = 2
                       , folds = 5
                       , spatial_folds = TRUE
                       , repeats = 1
                       , folds_div = seq(2.1, by = 0.1, length.out = folds)
                       , max_repeat_corr = 0.9
                       , min_fold_n = 8
                       , hold_prop = 0.3
                       , stretch_value = 10
                       , dens_res = 1000
                       , reduce_env_thresh_corr = 0.9
                       , reduce_env_quant_rf_imp = 0.2
                       , do_gc = FALSE
                       , force_new = FALSE
                       ) {

    # start presences -------
    start_presences <- if(pres_col %in% names(presence)) {

      sum(presence[pres_col] == pres_val)

    } else nrow(presence)


    # setup -------
    return_val <- if(any(isFALSE(out_dir), return_val == "object")) "prep" else "prep_file"

    if(isFALSE(out_dir)) {

      out_dir <- tempfile()

      delete_out <- TRUE

    } else {

      delete_out <- FALSE

    }

    ## out_dir ------
    # everything should go through here, at least as far as blah_file.rds
    if(is.character(out_dir)) {

      fs::dir_create(out_dir)

      if(dir.exists(out_dir)) {

        prep_file <- fs::path(out_dir
                              , "prep.rds"
                              )

        if(all(file.exists(prep_file), ! force_new)) {

          safe_import <- purrr::safely(rio::import)

          prep <- safe_import(prep_file, trust = TRUE)

          if(is.null(prep$error)) prep <- prep$result else {

            # remove the prep file if it can't be opened
            fs::file_delete(prep_file)

            rm(prep)

          }

        }

      } else stop("can't create out_dir")

    }

    if(!exists("prep", inherits = FALSE)) prep <- list(abandoned = FALSE, finished = FALSE, log = NULL)

    # run if not finished or abandoned, otherwise only run if force_new
    run <- if(any(prep$abandoned, prep$finished)) force_new else TRUE

    # but still don't run if start_presences < min_fold_n
    run <- if(start_presences < min_fold_n) FALSE else run

    if(is.null(this_taxa)) this_taxa <- basename(out_dir)

    prep$this_taxa <- this_taxa

    m <- paste0("prep for "
                , this_taxa
                , "\nout_dir is "
                , out_dir
                , ".\n "
                , start_presences
                , " incoming presences"
                )

    message(m)

    prep$log <- paste0(prep$log
                       , "\n"
                       , m
                       )

    if(run) {

      ## start timer -----
      start_time <- Sys.time()

      prep$log <- paste0(prep$log
                         , "\n"
                         , this_taxa
                         , "\nprep start at "
                         , start_time
                         , ".\n"
                         , nrow(presence)
                         , " rows in 'presence'"
                         , if(pres_col %in% names(presence)) paste0("\n of which "
                                                                    , presence |>
                                                                      dplyr::filter(!!rlang::ensym(pres_col) == pres_val) |>
                                                                      nrow()
                                                                    , " are presences."
                                                                    )
                         )

      ## predictors -----
      if(is_env_pred) {

        pred_names <- envRaster::name_env_tif(tibble::tibble(path = predictors), parse = TRUE) |>
          dplyr::pull(name)

        prep_preds <- terra::rast(predictors)
        names(prep_preds) <- pred_names

      } else {

        prep_preds <- terra::rast(predictors)

      }

      prep$log <- paste0(prep$log
                         , "\n"
                         , length(names(prep_preds))
                         , " predictors"
                         )

      # presence --------
      ## original -------
      prep$original <- presence

      ## raster pa -------
      prep$pa_ras <- terra::cellFromXY(prep_preds
                                       , presence %>%
                                         sf::st_as_sf(coords = c(pres_x, pres_y)
                                                      , crs = pres_crs
                                                      ) %>%
                                         sf::st_transform(crs = sf::st_crs(prep_preds)) %>%
                                         sf::st_coordinates()
                                       ) %>%
        terra::xyFromCell(prep_preds, .) %>%
        tibble::as_tibble() %>%
        {if(pres_col %in% names(presence)) (.) %>% dplyr::mutate(pa = presence$pa) else (.) |> dplyr::mutate(!!rlang::ensym(pres_col) := pres_val)} %>%
        stats::na.omit()

      ## raster presence ------
      # not limited to raster yet though
      prep$presence_ras <- prep$pa_ras |>
        dplyr::filter(!!rlang::ensym(pres_col) == pres_val) |>
        dplyr::select(- !!rlang::ensym(pres_col))


      # predict_boundary -------

      ## TRUE ------
      # create new mcp from points
      if(isTRUE(pred_limit)) {

        prep$predict_boundary <- prep$presence_ras %>%
          dplyr::rename(!!rlang::ensym(pres_x) := x
                        , !!rlang::ensym(pres_y) := y
                        ) |>
          dplyr::distinct(!!rlang::ensym(pres_x), !!rlang::ensym(pres_y)) %>%
          sf::st_as_sf(coords = c(pres_x, pres_y)
                       , crs = pres_crs
                       ) %>%
          sf::st_transform(crs = sf::st_crs(prep_preds[[1]])) %>%
          sf::st_union() %>%
          sf::st_convex_hull() %>%
          sf::st_sf() %>%
          sf::st_buffer(limit_buffer) %>%
          sf::st_make_valid()

      }

      ## file path -------
      # use existing MCP polygon file for the predict boundary
      if(is.character(pred_limit) & length(pred_limit)) {

        if(file.exists(pred_limit)) {

          prep$predict_boundary <- sfarrow::st_read_parquet(pred_limit) %>%
            sf::st_geometry() %>%
            sf::st_sf() %>%
            terra::vect() |>
            terra::densify(50000) |>
            sf::st_as_sf() |>
            sf::st_transform(crs = sf::st_crs(prep_preds[[1]])) %>%
            sf::st_make_valid()

        }

      }

      ## sf --------
      # existing sf
      if("sf" %in% class(pred_limit)) {

        prep$predict_boundary <- pred_limit %>%
          terra::vect() |>
          terra::densify(50000) |>
          sf::st_as_sf() |>
          sf::st_transform(crs = sf::st_crs(prep_preds[[1]])) %>%
          sf::st_make_valid() %>%
          sf::st_intersection(sf::st_bbox(prep_preds) %>%
                                sf::st_as_sfc() %>%
                                sf::st_sf()
                              ) %>%
          sf::st_make_valid()

      }

      # recast pred_limit from here on as T/F
      pred_limit <- if(all(!isFALSE(pred_limit), length(pred_limit))) TRUE else FALSE

      ## FALSE ------
      # use the full extent of the predictors for the predict boundary
      # No point buffering here
      # captures anything not already captured above
      if(isFALSE(pred_limit)) {

        prep$predict_boundary <- sf::st_bbox(prep_preds) %>%
          sf::st_as_sfc() %>%
          sf::st_sf()  %>%
          sf::st_make_valid()

      }

      # clip predict_boundary? -------
      if(!is.null(pred_clip)) {

        prep$predict_boundary <- prep$predict_boundary %>%
          sf::st_intersection(pred_clip %>%
                                sf::st_transform(crs = sf::st_crs(prep$predict_boundary)) |>
                                sf::st_make_valid()
                              ) %>%
          sf::st_make_valid()

      }

      # predict_boundary only on env rasters
      prep$predict_boundary <- prep$predict_boundary %>%
        sf::st_intersection(sf::st_bbox(prep_preds) %>%
                              sf::st_as_sfc() %>%
                              sf::st_sf()
                            )


      # prep$presence_ras clip to predict_boundary ---------
      # catch some cases where there are presence records on the predictors
      # but outside the predict boundary.
      prep$presence_ras <- prep$presence_ras %>%
        sf::st_as_sf(coords = c("x", "y")
                     , crs = sf::st_crs(prep_preds)
                     ) %>%
        sf::st_filter(prep$predict_boundary |>
                        terra::vect() |>
                        terra::densify(50000) |>
                        sf::st_as_sf() |>
                        sf::st_transform(crs = sf::st_crs(prep_preds)) |>
                        sf::st_make_valid()
                      )

      if(nrow(prep$presence_ras)) {

        prep$presence_ras <- prep$presence_ras %>%
          sf::st_coordinates() %>%
          tibble::as_tibble() %>%
          purrr::set_names(c("x", "y"))

      }

      # n_p / needed_p ---------
      n_p <- nrow(prep$presence_ras)
      needed_p <- min_fold_n * ((hold_prop > 0) + 1)

      prep$log <- paste0(prep$log
                         , "\n"
                         , needed_p
                         , " presences required for a single fold of "
                         , min_fold_n
                         , if(hold_prop > 0) paste0(" plus another fold of "
                                                    , min_fold_n
                                                    , " for holdout data"
                                                    )
                         , ". "
                         , n_p
                         , " presences in predict_boundary and on env rasters"
                         )

      if(n_p >= needed_p) {

        # subset predictors ------
        if(pred_limit) {

          if(!is.null(terra_options)) {

            do.call(terra::terraOptions
                    , args = terra_options
                    )

          }

          terra::window(prep_preds) <- terra::ext(terra::vect(prep$predict_boundary))

        }

        # folds adj -------

        # aiming for at least 'min_fold_n' presences in each fold

        k_folds <- folds

        while(ceiling(n_p * (1 - hold_prop)) < min_fold_n * k_folds) {

          k_folds <- k_folds - 1

        }

        if(k_folds < folds) {

          prep$log <- paste0(prep$log
                             , "\n"
                             , "warning: too few records to support original folds ("
                             , folds
                             , "). Folds reduced to "
                             , k_folds
                             )
        }

        k_folds_adj <- k_folds # use k_folds_adj if need to revert to non-spatial cv later (see folds/non-spatial)

        # num_bg adj------

        if(prop_abs == "prop") {

          num_bg <- num_bg * n_p

        }

        if(num_bg < many_p_prop * n_p) num_bg <- many_p_prop * n_p

        # density raster ------

        density_file <- fs::path(out_dir
                                 , "density.tif"
                                 )

        run <- if(file.exists(density_file)) force_new else TRUE

        if(run) {

          start_dens_ras <- Sys.time()

          # limit prep$pa_ras to predict_boundary
          prep$pa_ras <- prep$pa_ras |>
            envClean::filter_geo_range(use_aoi = prep$predict_boundary
                                       , x = "x"
                                       , y = "y"
                                       , crs_df = sf::st_crs(prep$predict_boundary)
                                       )

          if(all(!is.null(dens_res), !terra::is.lonlat(prep_preds[[1]]))) {

            # resolution of density raster < pred raster

            use_res <- if(terra::res(prep_preds)[[1]] < dens_res) dens_res else terra::res(prep_preds)[[1]]

            temp_ras <- terra::rast(resolution = use_res
                                    , crs = terra::crs(prep_preds[[1]])
                                    , extent = terra::ext(prep$predict_boundary)
                                    , vals = 1
                                    )

            bw <- MASS::kde2d(as.matrix(prep$pa_ras[,1])
                              , as.matrix(prep$pa_ras[,2])
                              , n = c(nrow(temp_ras), ncol(temp_ras))
                              , lims = terra::ext(temp_ras) %>% as.vector()
                              )

            target_density <- raster::raster(bw) %>%
              terra::rast() %>%
              terra::resample(temp_ras) %>%
              terra::focal(3
                           , mean
                           , na.policy = "only"
                           , na.rm = TRUE
                           )

          } else {

            # density raster = predict raster
            # STU METHOD
            # weighted bandwidth
            ## same method as spatialEco::sf.kde

            temp_ras <- terra::rast(resolution = terra::res(prep_preds)
                                    , crs = terra::crs(prep_preds[[1]])
                                    , extent = terra::ext(prep$predict_boundary)
                                    , vals = 1
                                    )

            pres_ras <- terra::rasterize(prep$pa_ras %>%
                                           sf::st_as_sf(coords = c("x", "y")
                                                        , crs = sf::st_crs(prep_preds[[1]])
                                                        )
                                         , y = temp_ras
                                         , fun = length
                                         , touches = TRUE
                                         )

            pres <- terra::as.data.frame(pres_ras, xy = TRUE) %>%
              tibble::as_tibble()

            colnames(pres)[3] <- "COUNT"

            pres <- pres %>%
              dplyr::mutate(scaled = COUNT * (length(COUNT) / sum(COUNT)))

            if(all(pres$scaled == 1)) pres$scaled <- sample(c(0.999, 1.001), length(pres$scaled), replace = TRUE) # nw

            bw <- ks::Hpi.diag(x = as.matrix(pres[, c("x", "y", "scaled")])
                               , pilot = "dscalar"
                               )

            target_density <- terra::setValues(pres_ras
                                               , matrix(ks::kde(as.matrix(pres[, c("x", "y")])
                                                                , eval.points = terra::xyFromCell(pres_ras, 1:terra::ncell(pres_ras))
                                                                , gridsize = c(nrow(pres_ras), ncol(pres_ras))
                                                                # use bw or NULL if error
                                                                ## NULL will use default method - maybe unweighted?
                                                                , h = ifelse(inherits(bw, "error"), NULL, bw)
                                                                , w = pres[["scaled"]]
                                                                , density = TRUE
                                                                )$estimate
                                                        , nrow = nrow(pres_ras)
                                                        , ncol = ncol(pres_ras)
                                                        , byrow = TRUE
                                                        )
                                               ) %>%
              terra::focal(3
                           , mean
                           , na.policy = "only"
                           , na.rm = TRUE
                           )

          }

          if(pred_limit) {

            target_density <- target_density %>%
              terra::mask(terra::vect(prep$predict_boundary)
                            , touches = TRUE
                          )

          }

          ## stretch-------
          target_density <- target_density %>%
            terra::stretch(1, stretch_value) %>%
            terra::subst(NA, 0)

          if(FALSE) {

            terra::plot(prep_preds[[1]])
            terra::plot(target_density
                        , add = TRUE
                        )

            ps <- prep$presence_ras %>%
              tibble::as_tibble() %>%
              sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(prep_preds)) %>%
              terra::vect()

            terra::plot(ps, add = TRUE)

            terra::plot(prep$predict_boundary %>% terra::vect()
                        , add = TRUE
                        )

            terra::plot(prep$bg_points %>% terra::vect()
                        , add = TRUE
                        )

          }

          terra::writeRaster(target_density
                             , density_file
                             , overwrite = TRUE
                             )

          prep$log <- paste0(prep$log
                             , "\n"
                             , "density raster done in "
                             , round(difftime(Sys.time(), start_dens_ras, units = "mins"), 2)
                             , " minutes"
                             )

        }


        # background points -------
        start_bg <- Sys.time()

        run <- if(exists("bg_points", where = prep)) force_new else TRUE

        if(run) {

          if(!exists("target_density")) target_density <- terra::rast(density_file)

          ptscell <- sample(1:terra::ncell(target_density)
                           , num_bg * 2 # over generate here as terra::window only crops - it does not mask
                           , prob = target_density[]
                           , replace = TRUE
                           )

          centres <- terra::xyFromCell(target_density, ptscell)

          hs <- terra::res(target_density) / 2

          prep$bg_points <- cbind(runif(nrow(centres), centres[, 1] - hs[1], centres[, 1] + hs[1])
                                  , runif(nrow(centres), centres[, 2] - hs[2], centres[, 2] + hs[2])
                                  ) %>%
            as.matrix() %>%
            terra::cellFromXY(prep_preds[[1]], .) %>%
            terra::xyFromCell(prep_preds[[1]], .) %>%
            tibble::as_tibble() %>%
            na.omit() %>%
            dplyr::distinct() %>%
            dplyr::anti_join(prep$presence_ras) %>% # background not on presences
            sf::st_as_sf(coords = c("x", "y")
                         , crs = sf::st_crs(target_density)
                         , remove = FALSE
                         ) %>%
            # filter back to only those within predict_boundary
            sf::st_filter(prep$predict_boundary)

          if(nrow(prep$bg_points) > num_bg) {

            prep$bg_points <- prep$bg_points %>%
              dplyr::sample_n(num_bg)

          }

          prep$log <- paste0(prep$log
                             , "\n"
                             , num_bg
                             , " background points attempted. Completed in: "
                             , round(difftime(Sys.time(), start_bg, units = "mins"), 2)
                             , " minutes"
                             )

        } else {

          target_density <- terra::rast(density_file)

        }

        # env--------
        start_env <- Sys.time()

        run <- if(exists("folds", prep)) force_new else TRUE

        if(run) {

          spp_pa <- dplyr::bind_rows(prep$presence_ras %>%
                                       tibble::as_tibble() %>%
                                       sf::st_as_sf(coords = c("x", "y")
                                                    , crs = sf::st_crs(prep_preds[[1]])
                                                    , remove = FALSE
                                                    ) %>%
                                       dplyr::mutate(pa = 1)
                                     , prep$bg_points %>%
                                       dplyr::mutate(pa = 0)
                                     )

          prep$env <- terra::extract(prep_preds
                                     , y = terra::vect(spp_pa)
                                     , include_cols = "pa"
                                     , ID = FALSE
                                     , bind = TRUE
                                     ) %>%
            tibble::as_tibble() |>
            na.omit()

          if(!is.null(cat_preds)) {

            prep$env <- prep$env %>%
              dplyr::mutate(dplyr::across(tidyselect::any_of(cat_preds), as.factor))

          }

          prep$log <- paste0(prep$log
                             , "\n"
                             , "env data extracted in: "
                             , round(difftime(Sys.time(), start_bg, units = "mins"), 2)
                             , " minutes.\n"
                             , nrow(prep$env[prep$env$pa != 1,])
                             , " background points with env values."
                             )

        }

        # Abandon if too few presences with env data
        if(nrow(prep$env[prep$env$pa == 1,]) < needed_p) {

          prep$log <- paste0(prep$log
                             , "\n"
                             , "warning: too few presences ("
                             , nrow(prep$env[prep$env$pa == 1,])
                             , ") with environmental variables. SDM abandoned"
                             )

          prep$abandoned <- TRUE

        }

        # folds------

        if(!prep$abandoned) {

          start_folds <- Sys.time()

          run <- if(exists("folds", prep)) force_new else TRUE

          if(run) {

            # testing/training split ----------
            to_split <- prep$env |>
              dplyr::mutate(id = dplyr::row_number())

            ## repeats adj -----
            repeats_adj <- if(n_p / min_fold_n < repeats) floor(n_p / min_fold_n) else repeats

            prep$testing <- tibble::tibble(pa = 0)
            counter <- 0
            hold_prop_adj <- hold_prop
            if(hold_prop_adj >= 0.5) hold_prop_adj = 0.49

            # Try to get min_fold_n presences within prep$testing
            if(hold_prop > 0) {

              while(any(sum(prep$testing$pa == 1) < min_fold_n, hold_prop_adj >= 0.5)) {

                hold_prop_adj <- hold_prop_adj + counter * 0.01

                prep$testing <- to_split %>%
                  dplyr::slice_sample(prop = hold_prop_adj, by = pa)

                counter <- counter + 1

              }

            } else {

              prep$testing <- to_split

            }

            prep$testing <- prep$testing |>
              dplyr::mutate(hold_prop = hold_prop_adj)

            if(hold_prop_adj > hold_prop) {

              m <- paste0("Warning. hold_prop adjusted up to "
                          , hold_prop_adj
                          , " (from "
                          , hold_prop
                          , ") to acheive min_fold_n of "
                          , min_fold_n
                          , " presences in the testing (holdout) data"
                          )

              message(m)

              prep$log <- paste0(prep$log
                                 , "\n"
                                 , m
                                 )

            }

            # If never enough presences in test - abandon
            if(sum(prep$testing$pa == 1) < min_fold_n) {

              m <- paste0("ERROR: not enough presences ("
                          , n_p
                          , ") to acheive min_fold_n ("
                          , min_fold_n
                          , ") in the test split with hold_prop of "
                          , hold_prop_adj
                          , " (adjusted up from "
                          , hold_prop
                          , ")"
                          )

              message(m)

              prep$log <- paste0(prep$log
                                 , "\n"
                                 , m
                                 )

              prep$abandoned <- TRUE

            }

          }

        }

        if(!prep$abandoned) {

          to_split <- to_split |>
            dplyr::anti_join(prep$testing |>
                               dplyr::distinct(id)
                             )

          if(! as.logical(nrow(to_split))) to_split <- prep$testing

          prep$training <- tibble::tibble(rep = 1:repeats_adj
                                          , training = list(to_split)
                                          ) |>
            dplyr::mutate(spatial_folds = dplyr::if_else(k_folds == 1, FALSE, spatial_folds))

          rm(to_split)

          text_one <- paste0("test/training split\n "
                             , as.character(nrow(prep$testing))
                             , " test data, including "
                             , as.character(sum(prep$testing$pa == 1))
                             , " presences."
                             )

          text_two <- paste0(" repeat: "
                             , 1:repeats_adj
                             , ". "
                             , purrr::map_chr(prep$training$training, \(x) as.character(nrow(x)))
                             , " training data, including "
                             , purrr::map_chr(prep$training$training, \(x) as.character(sum(x$pa == 1)))
                             , " presences."
                             )

          text <- paste0(c(text_one, text_two), collapse = "\n")

          prep$log <- paste0(prep$log
                             , "\n"
                             , text
                             )

          if(any(hold_prop == 0
                 , if(nrow(prep$training) == 1) prep$training$single_fold[[1]] else FALSE
                 )
             ) {

            prep$log <- paste0(prep$log
                               , "\n"
                               , "WARNING: no holdout data\n"
                               , " full model will be tested against the same data that was used to train it!"
                               )

          }

          if(k_folds > 1) {

          ## spatial ------
          safe_cv_spatial <- purrr::safely(blockCV::cv_spatial)

            prep$training <- prep$training |>
              dplyr::mutate(point_sf = purrr::map(training
                                                  , \(x) x |>
                                                    dplyr::select(x, y, pa) %>%
                                                    sf::st_as_sf(coords = c("x", "y")
                                                                 , crs = sf::st_crs(prep_preds[[1]])
                                                                 )
                                                  )
                            , folds_div = folds_div[1:nrow(prep$training)]
                            , fold_dist = purrr::map_dbl(folds_div
                                                          , \(x) prep$predict_boundary %>%
                                                            sf::st_area() %>%
                                                            as.numeric() %>%
                                                            sqrt() %>%
                                                            `/` (x)
                                                          )
                            , cv_spatial = purrr::map2(point_sf
                                                       , fold_dist
                                                       , \(x, y) suppressWarnings(
                                                         # The warnings are dealt with below by 'fix_folds'
                                                         safe_cv_spatial(x
                                                                         , column = "pa"
                                                                         , k = k_folds
                                                                         , size = y
                                                                         , iteration = 200
                                                                         , selection = "random"
                                                                         , extend = 0.5
                                                                         , progress = FALSE
                                                                         , report = FALSE
                                                                         , plot = FALSE
                                                                         )
                                                         )
                                                       )
                            , error = purrr::map(cv_spatial, "error")
                            )

            ### spatial cv errors -------
            if(any(purrr::map_lgl(prep$training$error, \(x) !is.null(x)))) {

              errs <- prep$training |>
                dplyr::filter(purrr::map_lgl(error, \(x) !is.null(x)))

              prep$log <- paste0(prep$log
                                 , "\n"
                                 , "error: repeat"
                                 , errs$rep
                                 , ". "
                                 , as.character(errs$error)
                                 )

            }

            ### extract folds from cv_spatial --------
            prep$training <- prep$training |>
              dplyr::filter(purrr::map_lgl(error, \(x) is.null(x))) |> # remove errors
              dplyr::mutate(folds = purrr::map(cv_spatial, \(x) x$result$folds_ids))

            ### spatial cv correlation --------
            if(nrow(prep$training) > 1) {

              reps <- prep$training$rep

              check_pres_corr <- purrr::map2(prep$training$folds
                                             , prep$training$training
                                             , \(x, y) x[y$pa == 1]
                                             ) |>
                purrr::set_names(reps)

              # catch any folds coming out of blockCV::cv_spatial that have all the presences in just one fold
              reps_single_fold <- purrr::map_lgl(check_pres_corr
                                                  , \(x) length(table(x)) == 1
                                                  ) |>
                purrr::set_names(reps)

              # run correlation between repeats (if there are at least 2 repeats with presences split among folds/folds)
              if(length(reps) - sum(reps_single_fold) > 1) {

                prep_fold_corr <- stats::cor(tibble::as_tibble(check_pres_corr[! reps_single_fold]
                                                               , .name_repair = "unique_quiet"
                                                               )
                                             )

                fold_corr_drop <- caret::findCorrelation(prep_fold_corr
                                                          , cutoff = max_repeat_corr
                                                          , names = TRUE
                                                          )

              } else fold_corr_drop <- NULL

              # combine single fold repeats and correlated repeats
              change_to_non_spatial <- c(names(reps_single_fold[reps_single_fold])
                                         , fold_corr_drop
                                         ) |>
                unique() |>
                as.numeric() |>
                sort()

              if(any(as.logical(change_to_non_spatial))) {

                prep$training$spatial_folds[change_to_non_spatial] <- FALSE

              }

            }

          }

          ## non-spatial option -------
          prep$training <- prep$training %>%
            dplyr::mutate(nsb = purrr::map(training
                                           , \(x) non_spatial_folds(k_folds
                                                                     , x
                                                                     )
                                           )
                          , folds = if(!"folds" %in% names(.)) {
                            nsb
                          } else {
                            folds
                          }
                          , folds = purrr::pmap(list(spatial_folds
                                                     , folds
                                                     , nsb
                                                     )
                                                , \(x, y, z) {

                                                  if(x) y else z

                                                }
                                                )
                          )

          ### too few p --------
          prep$training <- prep$training |>
            dplyr::mutate(folds = purrr::map2(folds, training
                                                 , \(x, y) fix_folds(x, y$pa, min_fold_n)
                                                 )
                          , single_fold = purrr::map_lgl(folds
                                                         , \(x) length(unique(x)) == 1
                                                         )
                          # again, if only one fold left, go to non-spatial
                          , folds = purrr::pmap(list(single_fold
                                                     , folds
                                                     , nsb
                                                     )
                                                , \(x, y, z) {

                                                  if(!x) y else z

                                                }
                                                )
                          )

          prep$training <- prep$training |>
            dplyr::mutate(training = purrr::map2(training
                                                 , folds
                                                 , \(x, y) x |>
                                                   dplyr::bind_cols(tibble::tibble(fold = y))
                                                 )
                          )

          # Find correlation for log
          if(length(prep$training$folds) > 1) {

            check_pres_corr <- purrr::map2(prep$training$folds
                                           , prep$training$training
                                           , \(x, y) x[y$pa == 1]
                                           )

            prep$prep_fold_corr <- stats::cor(tibble::as_tibble(check_pres_corr, .name_repair = "unique_quiet"))

            fold_corr <- prep$prep_fold_corr

            diag(fold_corr) <- NA

            max_corr <- max(fold_corr, na.rm = TRUE)

          }

          prep$log <- paste0(prep$log
                             , "\n"
                             , paste0("repeat: "
                                      , prep$training$rep
                                      , ", spatial folds = "
                                      , prep$training$spatial_folds
                                      , ". Folds = "
                                      , purrr::map(prep$training$folds
                                                   , \(x) x |>
                                                     unique() |>
                                                     length() |>
                                                     stringr::str_flatten_comma()
                                                   )
                                      , " with n presences = "
                                      , purrr::map(prep$training$training
                                                   , \(x) x |>
                                                     dplyr::count(fold, pa) |>
                                                     dplyr::filter(pa == 1) |>
                                                     dplyr::pull(n) |>
                                                     stringr::str_flatten_comma()
                                                   )
                                      , collapse = "\n"
                                      )
                             , if(length(prep$training$folds) > 1) paste0("\n Maximum correlation between repeats (presences only): "
                                                                           , round(max_corr, 3)
                                                                           )
                             , ".\nFolds done in "
                             , round(difftime(Sys.time(), start_folds, units = "mins"), 2)
                             , " minutes"
                             )

          # reduce env -------

          if(!prep$abandoned) {

            run <- if(exists("reduce_env", prep)) force_new else TRUE

            if(run) {

              start_reduce_env <- Sys.time()

              prep$reduce_env <- envModel::reduce_env(env_df = prep$env
                                                      , env_cols = names(prep_preds)
                                                      , y_col = "pa"
                                                      , imp_col = "1"
                                                      , thresh_corr = reduce_env_thresh_corr
                                                      , quant_rf_imp = reduce_env_quant_rf_imp
                                                      , remove_always = c(pres_x, pres_y, "x", "y", "pa", "fold", "cell", "id")
                                                      )

              prep$log <- paste0(prep$log
                                 , "\n"
                                 , "reduce_env completed in: "
                                 , round(difftime(Sys.time(), start_reduce_env, units = "mins"), 2)
                                 , " minutes. Original "
                                 , length(names(prep_preds))
                                 , " variables reduced to "
                                 , length(prep$reduce_env$keep)
                                 )

            }

          }

          # end timer ------

          prep$log <- paste0(prep$log
                             , "\n"
                             , "prep completed. elapsed time: "
                             , round(difftime(Sys.time(), start_time, units = "mins"), 2)
                             , " minutes"
                             )

          prep$finished <- TRUE

        }

      } else {

        prep$abandoned <- TRUE

        prep$log <- paste0(prep$log
                           , "\n"
                           , "only "
                           , n_p
                           , " useable presence points. SDM abandoned"
                           )

      }

      if(delete_out) {

        fs::dir_delete(out_dir)

      } else {

        rio::export(prep, prep_file)

      }

      # clean up -------
      if(do_gc) {

        stuff <- ls()

        stuff <- stuff[! stuff %in% c(return_val, "return_val")]

        rm(list = stuff)
        gc()

      }

    } else {

      # Should only get here if:
        # * prep.rds exists & (prep$abandoned or prep$finished) & force_new = FALSE
        # * there were less than min_fold_n presences

      if(start_presences < min_fold_n) {

        m <- paste0("Only "
                    , start_presences
                    , " incoming presences. SDM abandoned"
                    )

        message(m)

        prep$log <- paste0(prep$log
                           , "\n"
                           , m
                           )

        prep$abandoned <- TRUE
        prep$finished <- TRUE

        rio::export(prep, prep_file)

      }

    }

    # return-------
    res <- if(return_val == "prep") prep else prep_file

    return(res)

  }

