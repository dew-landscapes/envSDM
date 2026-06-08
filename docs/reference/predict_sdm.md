# Predict from SDM

The resulting `pred.tif` is masked to the boundary provided to the
`pred_limit` argument of prep_sdm; or generated in prep_sdm from the
`pred_limit`, `limit_buffer` and `pred_clip` arguments.

## Usage

``` r
predict_sdm(
  prep,
  full_run,
  out_dir,
  file_name = "pred.tif",
  use_env_naming = FALSE,
  predictors = NULL,
  is_env_pred = FALSE,
  terra_options = NULL,
  doClamp = TRUE,
  force_new = FALSE,
  do_gc = FALSE,
  check_tifs = FALSE,
  handle_errors = TRUE,
  ...
)
```

## Arguments

- prep:

  Character or named list. If character, the path to an existing
  `prep.rds`. Otherwise, the result of a call to
  [`prep_sdm()`](https://acanthiza.github.io/envSDM/reference/prep_sdm.md)
  with return_val = "object".

- full_run:

  Character or named list. If character, the path to an existing
  `full_run.rds`. Otherwise, the result of a call to
  [`run_full_sdm()`](https://acanthiza.github.io/envSDM/reference/run_full_sdm.md)
  with return_val = "object".

- out_dir:

  Character. Name of directory into which `.tif`s will be saved. Will be
  created if it does not exist.

- file_name:

  Character. Name to give the output prediction .tif.

- use_env_naming:

  Logical. If `TRUE`, and `is_env_pred` is `TRUE`, naming will ignore
  `file_name` and instead generate a name matching `name_env_tif()` with
  `layer` being `this_taxa` from `prep` and `start_date` being the
  minimum available `start_date` from the predictors. `pred` appears
  between `this_taxa` and `start_date`.

- predictors:

  Character. Vector of paths to predictor `.tif` files.

- is_env_pred:

  Logical. Does the naming of the directory and files in `predictors`
  follow the pattern required by `envRaster::parse_env_tif()`?

- terra_options:

  Passed to
  [`terra::terraOptions()`](https://rspatial.github.io/terra/reference/terraOptions.html).
  e.g. list(memfrac = 0.6)

- doClamp:

  Passed to
  [`terra::predict()`](https://rspatial.github.io/terra/reference/predict.html)
  (which then passes as `...` to `fun`). Possibly orphaned from older
  envSDM?

- force_new:

  Logical. If output files already exist, should they be remade?

- do_gc:

  Logical. Run `base::rm(list = ls)` and
  [`base::gc()`](https://rdrr.io/r/base/gc.html) at end of function?
  Useful to keep RAM use down when running SDMs for many, many taxa,
  especially if done in parallel.

- check_tifs:

  Logical. Check if any output `.tif` files error on
  [`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html)
  and delete them if they do. Useful after a crash during predict.

- handle_errors:

  Logical. Use purrr::safely when predicting, enabling the capture of
  (m)any errors (which are then written to the log). Suggest turning off
  (i.e. `handle_errors = FALSE`) when running in a targets pipeline.

- ...:

  Passed to `...` in
  [`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html) -
  the last step in the `envSDM::predict_sdm` process. Used to provide
  additional arguments to
  [`terra::writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.html).

## Value

Character path to predicted file, usually 'pred.tif'. Output .tif and
.log, written to `out_dir`.

## Examples

``` r
out_dir <- file.path(system.file(package = "envSDM"), "examples")

# setup -------
data <- readRDS(fs::path(out_dir, "data.rds"))

# predictors -------
preds <- fs::dir_ls(fs::path(out_dir, "tif"))

# Best combo--------
## run full SDM --------
future::plan(future::multisession())

furrr::future_pwalk(list(data$prep
                  , data$tune
                  , data$out_dir
                  )
             , \(a, b, c) run_full_sdm(prep = a
                                       , tune = b
                                       , out_dir = c
                                       , use_metric = "combo"

                                       # passed to tune_sdm via dots
                                       , metrics_df = envSDM::sdm_metrics
                                       #, force_new = FALSE
                                       )
             )


## predict -------
furrr::future_pwalk(list(data$prep
                  , data$out_dir
                  )
             , \(a, b) predict_sdm(prep = a
                                   , full_run = fs::path(b, "full_run.rds")
                                   , out_dir = b
                                   , predictors = preds
                                   , check_tifs = TRUE
                                   #, force_new = FALSE
                                   )
             )
#> Warning: UNRELIABLE VALUE: Future (<unnamed-73>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-73> (1f85c0626a33ed47c69f761da56504cb-73); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-74>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-74> (1f85c0626a33ed47c69f761da56504cb-74); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-75>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-75> (1f85c0626a33ed47c69f761da56504cb-75); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-76>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-76> (1f85c0626a33ed47c69f761da56504cb-76); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-77>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-77> (1f85c0626a33ed47c69f761da56504cb-77); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-78>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-78> (1f85c0626a33ed47c69f761da56504cb-78); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-79>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-79> (1f85c0626a33ed47c69f761da56504cb-79); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-80>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-80> (1f85c0626a33ed47c69f761da56504cb-80); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-81>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-81> (1f85c0626a33ed47c69f761da56504cb-81); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-82>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-82> (1f85c0626a33ed47c69f761da56504cb-82); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-83>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-83> (1f85c0626a33ed47c69f761da56504cb-83); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-84>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-84> (1f85c0626a33ed47c69f761da56504cb-84); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-85>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-85> (1f85c0626a33ed47c69f761da56504cb-85); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-86>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-86> (1f85c0626a33ed47c69f761da56504cb-86); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-87>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-87> (1f85c0626a33ed47c69f761da56504cb-87); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-88>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-88> (1f85c0626a33ed47c69f761da56504cb-88); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-89>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-89> (1f85c0626a33ed47c69f761da56504cb-89); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-90>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-90> (1f85c0626a33ed47c69f761da56504cb-90); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-91>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-91> (1f85c0626a33ed47c69f761da56504cb-91); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-92>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-92> (1f85c0626a33ed47c69f761da56504cb-92); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-93>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-93> (1f85c0626a33ed47c69f761da56504cb-93); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-94>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-94> (1f85c0626a33ed47c69f761da56504cb-94); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-95>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-95> (1f85c0626a33ed47c69f761da56504cb-95); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-96>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-96> (1f85c0626a33ed47c69f761da56504cb-96); on 1f85c0626a33ed47c69f761da56504cb@bi-mod-01<2011444>]

future::plan(future::sequential())

## visualise-------
# just use one taxa
vis_data <- data |>
  dplyr::filter(taxa == "chg")

tifs <- fs::path(vis_data$out_dir[file.exists(fs::path(vis_data$out_dir, "pred.tif"))], "pred.tif")

names <- paste0("hold_prop "
                , vis_data$hold_prop
                , "; stretch "
                , vis_data$stretch
                , "; spatial_folds "
                , vis_data$spatial_folds
                )

r <- terra::rast(tifs)
names(r) <- names
terra::panel(r, cex.main = 0.6, nc = 2)
```
