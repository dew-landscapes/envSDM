
#' Prepare for running an SDM
#'
#' The background sampling includes code based on a
#' [Geographic Information Systems stack exchange](https://gis.stackexchange.com/)
#' [answer](https://gis.stackexchange.com/a/224347)
#' by user [Spacedman](https://gis.stackexchange.com/users/865/spacedman).
#'
#' @param this_taxa Character. Name of taxa. Used to name outputs. If `NULL`,
#' this will be `basename(dirname(out_dir))`.
#' @param out_dir Character. Name of directory into which results will be saved.
#' @param presence Dataframe of presences. Needs columns called `lat` (latitude)
#' and `long` (longitude) both in decimal degrees.
#' @param pres_crs Anything that will return a legitimate crs when passed to the
#' crs attribute of st_transform or st_as_sf
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates
#' @param pred_limit Limit the background points and predictions?
#' Can be `TRUE` (use points to generate a minimum convex polygon to use as a
#' limit), `FALSE` (do not limit) or path to existing sf to use.
#' @param limit_buffer Numeric. Apply this buffer to `pred_limit`. This value
#' is passed to the `dist` argument of `sf::st_buffer()`.
#' @param pred_clip sf. Optional sf to clip the pred_limit back to (e.g. to
#' prevent prediction into ocean)
#' @param predictors Character. Vector of paths to predictor `.tif` files.
#' @param is_env_pred Logical. Does the naming of the directory and files in
#' `predictors` follow the pattern required by `envRaster::parse_env_tif()`?
#' @param cat_preds Character. Vector of predictor names that are character.
#' @param num_bg Numeric. How many background points?
#' @param prop_abs Character. Is `num_bg` a proportion (`prop`) of the number of
#'  records in `presence` or an absolute (`abs`) number?
#' @param many_p_prop Numeric. Ensure the number of background points is at
#' least `many_p_prop * number of presences`. e.g. If there are more than 5000
#' presences and num_bg is set at `10000` and `many_p_prop` is `2`, then num_bg
#' will be increased to `many_p_prop * nrow(presences)`
#' @param folds Numeric. How many folds to use in cross validation? Will be
#' adjusted downwards if number of presences do not support `folds * min_fold_n`
#' @param spatial_folds Logical. Use spatial folds? Even if `TRUE`, can resort
#' to non-spatial cv if presences per fold do not meet `min_fold_n` or there are
#' not enough presences to support more than one fold.
#' @param min_fold_n Numeric. Sets both minimum number of presences, and,
#' by default, the minimum number of presences required for a model.
#' @param stretch_value Numeric. Stretch the density raster to this value.
#' @param dens_res `NULL` or numeric. Resolution (in metres) of density raster.
#' Set to `NULL` to use the same resolution as the predictors.
#' @param save_pngs Logical. Save out a .png of the density raster and spatial
#' blocks
#' @param remove_corr Logical. If TRUE, predictors with high correlation are
#' removed.
#' @param corr_thresh Numeric. Definition of 'high' correlation for `remove_corr`
#' @param do_gc Logical. Run `base::rm(list = ls)` and `base::gc()` at end of
#' function? Useful when running SDMs for many, many taxa, especially if done in
#' parallel.
#' @param force_new Logical. If outputs already exist, should they be remade?
#'
#' @return `invisible(NULL)`. `prep.rds` (a list) and log written to `out_dir`
#' @export
#'
#' @example inst/examples/prep_sdm_ex.R
  prep_sdm <- function(this_taxa = NULL
                       , out_dir
                       , presence
                       , pres_crs = 4326
                       , pres_x = "long"
                       , pres_y = "lat"
                       , pred_limit = TRUE
                       , limit_buffer = 0
                       , pred_clip = NULL
                       , predictors
                       , is_env_pred = TRUE
                       , cat_preds = NULL
                       , num_bg = 10000
                       , prop_abs = "abs"
                       , many_p_prop = 2
                       , folds = 5
                       , spatial_folds = TRUE
                       , min_fold_n = 8
                       , stretch_value = 10
                       , dens_res = 1000
                       , save_pngs = TRUE
                       , remove_corr = TRUE
                       , corr_thresh = 0.9
                       , do_gc = FALSE
                       , force_new = FALSE
                       ) {

    prep_log <- fs::path(out_dir
                           , "prep.log"
                           )

    run <- if(file.exists(prep_log)) {

      any(!grepl("prep end", paste0(readLines(prep_log), collapse = " "))
          , force_new
          )

    } else TRUE

    if(run) {

      if(is.null(this_taxa)) this_taxa <- basename(out_dir)

      message(paste0("prep for ", this_taxa))

      fs::dir_create(out_dir)

      # setup ------

      prep_file <- fs::path(out_dir
                            , "prep.rds"
                            )

      prep <- list(inputs = mget(ls(pattern = "[^predictors]"))
                   , original = presence
                   )

      ## start timer -----
      prep$timer <- envFunc::timer(process = "prep start"
                          , notes = paste0(nrow(presence), " original record"
                                           , if(nrow(presence) > 1) "s"
                                           )
                          , file = "prep"
                          , name = this_taxa
                          , log = prep_log
                          , write_log = TRUE
                          )

      ## predictors -----
      if(is_env_pred) {

        pred_names <- envRaster::name_env_tif(tibble::tibble(path = predictors), parse = TRUE) %>%
          dplyr::mutate(name = paste0(season,"__", layer)) %>%
          dplyr::pull(name)

        predictors <- terra::rast(predictors)
        names(predictors) <- pred_names

      } else {

        predictors <- terra::rast(predictors)

      }

      p <- presence %>%
        sf::st_as_sf(coords = c(pres_x, pres_y)
                     , crs = pres_crs
                     ) %>%
        sf::st_transform(crs = sf::st_crs(predictors)) %>%
        sf::st_coordinates()

      # raster presence ------

      prep$presence <- terra::cellFromXY(predictors
                                         , unique(p)
                                         ) %>%
        unique() %>%
        terra::xyFromCell(predictors, .) %>%
        stats::na.omit() %>%
        tibble::as_tibble()

      prep$timer <- envFunc::timer("raster presence"
                          , notes = paste0(nrow(prep$presence), " raster presences")
                          , time_df = prep$timer
                          , write_log = TRUE
                          )

      if(nrow(prep$presence) > min_fold_n) {

        # predict limits -------

        if(!isFALSE(pred_limit)) {

          if(is.character(pred_limit)) {
            # use existing MCP polygon file for the predict boundary

            prep$predict_boundary <- sfarrow::st_read_parquet(pred_limit) %>%
              sf::st_transform(crs = sf::st_crs(predictors[[1]])) %>%
              sf::st_buffer(limit_buffer) %>%
              sf::st_make_valid() %>%
              sf::st_sf()

          } else {
            # create new MCP polygon around the presences for the predict boundary

            prep$predict_boundary <- presence %>%
              dplyr::distinct(!!rlang::ensym(pres_x), !!rlang::ensym(pres_y)) %>%
              sf::st_as_sf(coords = c(pres_x, pres_y)
                           , crs = pres_crs
                           ) %>%
              sf::st_transform(crs = sf::st_crs(predictors[[1]])) %>%
              sf::st_union() %>%
              sf::st_convex_hull() %>%
              sf::st_sf() %>%
              sf::st_buffer(limit_buffer) %>%
              sf::st_make_valid()

          }

          pred_limit <- TRUE

        } else {
          # use the full extent of the predictors for the predict boundary

          prep$predict_boundary <- sf::st_bbox(predictors) %>%
            sf::st_as_sfc() %>%
            sf::st_sf()  %>%
            sf::st_make_valid()

        }

        if(!is.null(pred_clip)) {

          prep$predict_boundary <- prep$predict_boundary %>%
            sf::st_intersection(pred_clip %>%
                                  sf::st_transform(crs = sf::st_crs(prep$predict_boundary))
                                ) %>%
            sf::st_make_valid()


        }

        prep$timer <- envFunc::timer("predict boundary"
                            , time_df = prep$timer
                            , write_log = TRUE
                            )

        # prep$presence clip to predict_boundary ---------
        # catch some cases where there are presence records on the predictors
        # but outside the predict boundary (so there are no background points
        # around those presences).
        prep$presence <- prep$presence %>%
          sf::st_as_sf(coords = c("x", "y")) %>%
          sf::st_transform(crs = sf::st_crs(prep$predict_boundary)) %>%
          sf::st_intersection() %>%
          sf::st_set_geometry(NULL)


        # folds adj -------

        # aiming for at least 'min_fold_n' presences in each fold

        k_folds <- folds

        while(nrow(prep$presence) < min_fold_n * k_folds) {

          k_folds <- k_folds - 1

        }

        if(k_folds == 0) k_folds <- 1

        if(k_folds < folds) {

          prep$timer <- envFunc::timer(process = "warning"
                              , notes = paste0("too few records to support original folds (", folds, "). Folds reduced to ", k_folds)
                              , time_df = prep$timer
                              , write_log = TRUE
                              )

        }

        if(k_folds < 2) {

          prep$timer <- envFunc::timer(process = "warning"
                              , notes = paste0("WARNING: folds = ", k_folds)
                              , time_df = prep$timer
                              , write_log = TRUE
                              )

          spatial_folds <- FALSE

        }

        # num_bg adj------

        if(prop_abs == "prop") {

          num_bg <- num_bg * nrow(prep$presence)

        }

        if(num_bg < many_p_prop * nrow(prep$presence)) num_bg <- many_p_prop * nrow(prep$presence)

        # density raster ------

        density_file <- fs::path(out_dir
                                 , "density.tif"
                                 )

        run <- if(file.exists(density_file)) force_new else TRUE

        if(run) {

          if(all(!is.null(dens_res), !terra::is.lonlat(predictors[[1]]))) {

            # resolution of density raster < pred raster

            use_res <- if(terra::res(predictors)[[1]] < dens_res) dens_res else terra::res(predictors)[[1]]

            temp_ras <- terra::rast(resolution = use_res
                                    , crs = terra::crs(predictors[[1]])
                                    , extent = terra::ext(prep$predict_boundary)
                                    , vals = 1
                                    )

            bw <- MASS::kde2d(as.matrix(prep$presence[,1])
                              , as.matrix(prep$presence[,2])
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

            temp_ras <- terra::rast(resolution = terra::res(predictors)
                                    , crs = terra::crs(predictors[[1]])
                                    , extent = terra::ext(prep$predict_boundary)
                                    , vals = 1
                                    )

            pres_ras <- terra::rasterize(prep$presence %>%
                                           sf::st_as_sf(coords = c(x, y)
                                                        , crs = sf::st_crs(predictors[[1]])
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

            use_mask <- prep$predict_boundary %>%
              sf::st_transform(crs = sf::st_crs(target_density[[1]])) %>%
              sf::st_make_valid()

            target_density <- target_density %>%
              terra::mask(terra::vect(use_mask)
                            , touches = TRUE
                          )

          }

          ## stretch-------
          target_density <- target_density %>%
            terra::stretch(1, stretch_value) %>%
            terra::subst(NA,0)

          if(FALSE) terra::plot(target_density)

          terra::writeRaster(target_density
                             , density_file
                             , overwrite = TRUE
                             )

          prep$timer <- envFunc::timer("density"
                              , time_df = prep$timer
                              , write_log = TRUE
                              )

        }


        # background points -------
        run <- if(exists("bg_points", where = prep)) force_new else TRUE

        if(run) {

          if(!exists("target_density")) target_density <- terra::rast(density_file)

          ptscell <- sample(1:terra::ncell(target_density)
                           , num_bg * 1.1
                           , prob = target_density[]
                           , replace = TRUE
                           )

          centres <- terra::xyFromCell(target_density, ptscell)

          hs <- terra::res(target_density) / 2

          prep$bg_points <- cbind(runif(nrow(centres), centres[, 1] - hs[1], centres[, 1] + hs[1])
                                  , runif(nrow(centres), centres[, 2] - hs[2], centres[, 2] + hs[2])
                                  ) %>%
            as.matrix() %>%
            terra::cellFromXY(predictors[[1]], .) %>%
            terra::xyFromCell(predictors[[1]], .) %>%
            tibble::as_tibble() %>%
            na.omit() %>%
            dplyr::distinct() %>%
            sf::st_as_sf(coords = c("x", "y")
                         , crs = sf::st_crs(target_density)
                         )

          if(nrow(prep$bg_points) > num_bg) {

            prep$bg_points <- prep$bg_points %>%
              dplyr::sample_n(num_bg)

          }

          rio::export(prep, prep_file)

          prep$timer <- envFunc::timer("background points"
                                , notes = paste0(num_bg, " background points attempted. ", nrow(prep$bg_points), " achieved.")
                                , time_df = prep$timer
                                , write_log = TRUE
                                )


        } else {

        target_density <- terra::rast(density_file)

        prep <- rio::import(prep_file)

      }

      if(save_pngs) {

        # density png ------

        dens_png <- fs::path(out_dir, "density.png")

        if(!file.exists(dens_png)) {

          png_from_tif(target_density
                       , title = this_taxa
                       , dots = prep$presence %>%
                         sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(target_density))
                       , trim = TRUE
                       , out_png = dens_png
                       )

        }

      }


      # env--------

      run <- if(exists("spp_pa_env", prep)) force_new else TRUE

      if(run) {

        spp_pa <- dplyr::bind_rows(prep$presence %>%
                                     dplyr::select(x, y) %>%
                                     sf::st_as_sf(coords = c("x", "y")
                                                  , crs = sf::st_crs(predictors[[1]])
                                                  ) %>%
                                     dplyr::mutate(pa = 1)
                                   , prep$bg_points %>%
                                     dplyr::mutate(pa = 0) %>%
                                     dplyr::select(pa)
                                   ) %>%
          sf::st_buffer(terra::res(predictors)[[1]] / 100)

        prep$spp_pa_env <- exactextractr::exact_extract(predictors
                                                   , y = spp_pa
                                                   , include_cols = "pa"
                                                   , include_cell = TRUE
                                                   , include_xy = TRUE
                                                   , max_cells_in_memory = 40000000
                                                   ) %>%
          dplyr::bind_rows() %>%
          stats::na.omit() %>%
          dplyr::select(-coverage_fraction) %>%
          dplyr::distinct()

        if(!is.null(cat_preds)) {

          prep$spp_pa_env <- prep$spp_pa_env %>%
            dplyr::mutate(dplyr::across(tidyselect::any_of(cat_preds), as.factor))

        }

        rio::export(prep
                    , prep_file
                    )

        prep$timer <- envFunc::timer("env data"
                            , time_df = prep$timer
                            , write_log = TRUE
                            )

      }

      # Abandon if too few presences with env data
      if(nrow(prep$spp_pa_env[prep$spp_pa_env$pa == 1,]) < min_fold_n) {

        prep$timer <- envFunc::timer("warning"
                            , notes = paste0("Too few presences ("
                                             , nrow(prep$spp_pa_env[prep$spp_pa_env$pa == 1,])
                                             , ") with environmental variables. SDM abandoned"
                                             )
                            , time_df = prep$timer
                            )

        abandoned <- TRUE

      } else {

        abandoned <- FALSE

      }


      # folds------

      if(!abandoned) {

        run <- if(exists("blocks", prep)) force_new else TRUE

        if(run) {

          if(!all(spatial_folds, k_folds > 1)) {

            spatial_folds <- FALSE

          } else {

            ## spatial ------

            safe_cv_spatial <- purrr::safely(blockCV::cv_spatial)

            x <- prep$spp_pa_env %>%
              na.omit() %>%
              dplyr::select(x, y, pa) %>%
              sf::st_as_sf(coords = c("x", "y")
                           , crs = sf::st_crs(predictors[[1]])
                           )

            block_dist <- prep$predict_boundary %>%
              sf::st_as_sf() %>%
              dplyr::summarise() %>%
              sf::st_convex_hull() %>%
              sf::st_area() %>%
              as.numeric() %>%
              sqrt() %>%
              `/` (6)
            # for a square MCP, this would give 6 blocks by 6 blocks

            blocks <- safe_cv_spatial(x
                                      , column = "pa"
                                      #, r = predictors[[1]] # hashed out 2024-09-17
                                      , k = k_folds
                                      , size = block_dist
                                      , iteration = 200
                                      , selection = "random"
                                      , extend = 0.5
                                      , progress = FALSE
                                      , report = FALSE
                                      , plot = FALSE
                                      )

            if(!is.null(blocks$error)) {

              prep$timer <- envFunc::timer("error"
                                  , notes = paste0(as.character(blocks$error)
                                                   , ". spatial_folds set to FALSE"
                                                   )
                                  , time_df = prep$timer
                                  , write_log = TRUE
                                  )

              spatial_folds <- FALSE

              blocks <- NULL

            } else {

              blocks <- blocks$result$folds_ids %>%
                tibble::enframe(name = NULL, value = "fold_ids")

            }

            if(!is.null(blocks)) {

              blocks_p <- blocks[prep$spp_pa_env$pa == 1,]

              if(any(c(table(blocks_p$fold_ids) < min_fold_n), length(setdiff(1:k_folds, unique(blocks_p$fold_ids))) > 0)) {

                how_many_below_thresh <- sum(purrr::map_dbl(1:k_folds, ~ sum(blocks_p$fold_ids == .)) < min_fold_n)

                old_k_folds <- k_folds

                k_folds <- old_k_folds - how_many_below_thresh

                blocks_adj <- tibble::tibble(fold_ids = 1:old_k_folds) %>%
                  dplyr::mutate(n = purrr::map_dbl(fold_ids, ~ sum(blocks_p$fold_ids == .))) %>%
                  dplyr::mutate(fold_ids_adj = forcats::fct_lump_n(as.factor(fold_ids), k_folds - 1, w = n)) %>%
                  dplyr::mutate(fold_ids_adj = as.numeric(fold_ids_adj)) %>%
                  dplyr::distinct()

                blocks <- blocks %>%
                  dplyr::left_join(blocks_adj) %>%
                  dplyr::mutate(fold_ids = fold_ids_adj) %>%
                  dplyr::select(-fold_ids_adj, -n)

                k_folds <- length(unique(blocks$fold_ids))

                failed_blocks <- max(blocks_adj$fold_ids) - max(blocks_adj$fold_ids_adj)

                if(max(blocks_adj$fold_ids_adj) == 1) {

                  note <- paste0(failed_blocks
                                 , " out of "
                                 , nrow(blocks_adj)
                                 , " blocks failed to reach "
                                 , min_fold_n
                                 , " presences, leaving only 1 block. Reverting to non-spatial cv"
                                 )

                  spatial_folds <- FALSE

                } else {

                  note <- paste0(failed_blocks
                                 , " out of "
                                 , nrow(blocks_adj)
                                 , " blocks failed to reach "
                                 , min_fold_n
                                 , " presences. These were lumped until every block reached "
                                 , min_fold_n
                                 , " presences."
                                 )

                }

                prep$timer <- envFunc::timer("warning"
                                    , notes = note
                                    , time_df = prep$timer
                                    , write_log = TRUE
                                    )

              }

            }

          }

          if(any(!spatial_folds, k_folds == 1)) {

            ## non-spatial -------

            blocks <- tibble::tibble(fold_ids = c(sample(1:k_folds
                                                         , sum(prep$spp_pa_env$pa == 1)
                                                         , replace = TRUE
                                                         , prob = rep(1 / k_folds, k_folds)
                                                         )
                                                  , sample(1:k_folds
                                                           , sum(prep$spp_pa_env$pa == 0)
                                                           , replace = TRUE
                                                           , prob = rep(1 / k_folds, k_folds)
                                                           )
                                                  )
                                     )

          }

          prep$blocks <- prep$spp_pa_env %>%
            dplyr::mutate(block = blocks$fold_ids) %>%
            dplyr::filter(dplyr::if_any(.cols = names(predictors)
                                        , .fns = \(x) !is.na(x) & !is.infinite(x)
                                        )
                          )

          prep$spatial_folds_used <- spatial_folds

          prep$timer <- envFunc::timer("folds"
                              , notes = paste0("spatial folds = ", spatial_folds, ". Folds = ", length(unique(blocks$fold_ids)))
                              , time_df = prep$timer
                              , write_log = TRUE
                              )

          rio::export(prep
                      , prep_file
                      )


          ### block png -------
          if(save_pngs) {

            block_file <- fs::path(out_dir, "blocks.png")

            title <- paste0(this_taxa
                            , "\nSpatial folds = "
                            , prep$spatial_folds_used
                            , "\nFolds = "
                            , length(unique(prep$blocks$block))
                            , "\nPresences = "
                            , format(nrow(prep$blocks[prep$blocks$pa == 1,]), big.mark = ",")
                            , "\nBackground = "
                            , format(nrow(prep$blocks[prep$blocks$pa == 0,]), big.mark = ",")
                            )

            m <- prep$blocks %>%
              dplyr::select(x, y, block) %>%
              sf::st_as_sf(coords = c("x", "y")
                           , crs = sf::st_crs(prep$bg_points)
                           ) %>%
              tmap::tm_shape() +
                tmap::tm_dots(col = "block"
                              , size = 0.01
                              , palette = "viridis"
                              , breaks = seq(1, max(prep$blocks$block) + 1, 1)
                              , labels = as.character(seq(1, max(prep$blocks$block), 1))
                              , alpha = 0.5
                              ) +
              tmap::tm_credits(text = title
                               , position = c("left", "bottom")
                               )

              tmap::tmap_save(m
                              , block_file
                              )

          }

        }

      }


      # correlated -------

      if(all(remove_corr, !abandoned)) {

        run <- if(exists("correlated", prep)) force_new else TRUE

        if(run) {

          prep$correlated <- envModel::make_env_corr(prep$blocks
                                                     , env_cols = names(predictors)
                                                     , remove = TRUE
                                                     , thresh = corr_thresh
                                                     , always_remove = c(pres_x, pres_y, "x", "y", "pa", "block", "cell")
                                                     )

          prep$timer <- envFunc::timer("correlation"
                                       , time_df = prep$timer
                                       , write_log = TRUE
                                       )

        }

      } else {

        if(!abandoned) {

          prep$correlated$remove_env <- ""
          prep$correlated$env_cols <- names(predictors)

        }

      }

        # end timer ------
        prep$timer <- envFunc::timer(process = "prep end"
                                     , time_df = prep$timer
                                     )

        rio::export(prep
                  , prep_file
                  )

      } else {

        prep$timer <- envFunc::timer("warning"
                            , notes = paste0("only ", nrow(prep$presence), " (raster) presence points. SDM abandoned")
                            , time_df = prep$timer
                            )

        prep$timer <- envFunc::timer(process = "prep end"
                            , time_df = prep$timer
                            )

        rio::export(prep
                    , prep_file
                    )

      }

      if(do_gc) {

        rm(list = ls())

        gc()

      }

    }

    return(invisible(NULL))

  }

