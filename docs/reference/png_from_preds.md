# Generate .png (map) files for each prediction

Finds any .tif files in `pred_dir` and writes them to .png files.
Includes the retrieval and addition to the map of: various SDM metrics;
and the original presence points.

## Usage

``` r
png_from_preds(
  pred_dir,
  full_run_dir = NULL,
  trim = TRUE,
  force_new = FALSE,
  do_gc = TRUE,
  ...
)
```

## Arguments

- trim:

  Logical. Trim NA values from outside (using
  [`terra::trim()`](https://rspatial.github.io/terra/reference/trim.html))

- force_new:

  Logical. If .png file already exists, recreate it?

- do_gc:

  Logical. Run `base::rm(list = ls)` and
  [`base::gc()`](https://rdrr.io/r/base/gc.html) at end of function?
  Useful when running SDMs for many, many taxa, especially if done in
  parallel.

- ...:

  Passed to [`fs::dir_ls()`](https://fs.r-lib.org/reference/dir_ls.html)

- prep:

  Character or named list. If character, the path to an existing
  `prep.rds`. Otherwise, the result of a call to
  [`prep_sdm()`](https://acanthiza.github.io/envSDM/reference/prep_sdm.md)
  with return_val = "object".

- full_run:

  Character or named list. If character, the path to an existing
  `full_run.rds`. Otherwise, the result of a call to `run_full_sdm`()
  with return_val = "object".

- out_dir:

  Character. Name of directory into which `.pngs`s will be saved. Will
  be created if it does not exist.

## Value

`invisible(NULL)`. Writes .png files with the same file name as any .tif
files

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
#> Warning: UNRELIABLE VALUE: Future (<unnamed-25>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-25> (9220e353f90d2e3fdca4103a4d4c0f1f-25); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-26>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-26> (9220e353f90d2e3fdca4103a4d4c0f1f-26); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-27>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-27> (9220e353f90d2e3fdca4103a4d4c0f1f-27); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-28>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-28> (9220e353f90d2e3fdca4103a4d4c0f1f-28); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-29>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-29> (9220e353f90d2e3fdca4103a4d4c0f1f-29); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-30>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-30> (9220e353f90d2e3fdca4103a4d4c0f1f-30); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-31>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-31> (9220e353f90d2e3fdca4103a4d4c0f1f-31); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-32>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-32> (9220e353f90d2e3fdca4103a4d4c0f1f-32); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-33>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-33> (9220e353f90d2e3fdca4103a4d4c0f1f-33); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-34>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-34> (9220e353f90d2e3fdca4103a4d4c0f1f-34); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-35>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-35> (9220e353f90d2e3fdca4103a4d4c0f1f-35); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-36>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-36> (9220e353f90d2e3fdca4103a4d4c0f1f-36); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-37>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-37> (9220e353f90d2e3fdca4103a4d4c0f1f-37); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-38>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-38> (9220e353f90d2e3fdca4103a4d4c0f1f-38); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-39>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-39> (9220e353f90d2e3fdca4103a4d4c0f1f-39); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-40>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-40> (9220e353f90d2e3fdca4103a4d4c0f1f-40); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-41>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-41> (9220e353f90d2e3fdca4103a4d4c0f1f-41); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-42>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-42> (9220e353f90d2e3fdca4103a4d4c0f1f-42); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-43>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-43> (9220e353f90d2e3fdca4103a4d4c0f1f-43); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-44>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-44> (9220e353f90d2e3fdca4103a4d4c0f1f-44); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-45>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-45> (9220e353f90d2e3fdca4103a4d4c0f1f-45); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-46>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-46> (9220e353f90d2e3fdca4103a4d4c0f1f-46); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-47>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-47> (9220e353f90d2e3fdca4103a4d4c0f1f-47); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]
#> Warning: UNRELIABLE VALUE: Future (<unnamed-48>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-48> (9220e353f90d2e3fdca4103a4d4c0f1f-48); on 9220e353f90d2e3fdca4103a4d4c0f1f@bi-mod-01<800160>]

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
