# Create a .png from a .tif

Used within prep_sdm to save png of density raster with presences

## Usage

``` r
png_from_tif(
  x,
  title = NULL,
  dots = NULL,
  trim = TRUE,
  out_png = NULL,
  do_gc = FALSE
)
```

## Arguments

- x:

  spatRaster or path to .tif

- title:

  Character. Title to add to the .png

- dots:

  sf. Usually presences. Added as points.

- trim:

  Logical. Run
  [`terra::trim()`](https://rspatial.github.io/terra/reference/trim.html)
  before writing to .png?

- out_png:

  Character. Name of .png file to save. If `NULL` will be the same file
  name as `terra::sources(x)` with the file type as .png

- do_gc:

  Logical. Run `base::rm(list = ls)` and
  [`base::gc()`](https://rdrr.io/r/base/gc.html) at end of function?
  Useful when running SDMs for many, many taxa, especially if done in
  parallel.

## Value

`invisible(NULL)`. `out_png` written.
