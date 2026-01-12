# Merge polygons (or polygon files) to form a single minimum convex polygon (mcp)

Primary use is to create a predict boundary for
[`prep_sdm()`](https://acanthiza.github.io/envSDM/reference/prep_sdm.md),
merging (an existing) minimum convex polygon with other sources of taxa
distribution and clipping to a (usually coastal) boundary. The predict
boundary is then used for generation of background points and for
masking the 'full' predict (to the extent of the environmental
variables.

## Usage

``` r
make_predict_boundary(
  poly_list,
  buffer_list = NULL,
  out_file = tempfile(),
  col_name = "taxa",
  col_name_val = "boundary",
  clip = NULL,
  out_crs,
  return_poly = FALSE,
  force_new = FALSE,
  dens_int = NULL
)
```

## Arguments

- poly_list:

  List of paths or list of sf.

- buffer_list:

  List, of length 1 or the same length as `poly_list`, specifying the
  distance to buffer each incoming element of `poly_list`. Buffering
  occurs after transformation to `out_crs`.

- out_file:

  Character name of file to save

- col_name:

  Name of column to create in the resulting mcp

- col_name_val:

  Value to provide in the column in the resulting mcp

- clip:

  sf. Clip the resulting mcp back to this.

- out_crs:

  Numeric. [epsg](https://epsg.io/) code

- return_poly:

  Logical. Return the mcp, or alternatively `out_file`

- force_new:

  Logical. If `out_file` exists, recreate it?

- dens_int:

  Numeric.
  [`terra::densify()`](https://rspatial.github.io/terra/reference/densify.html)
  `interval` argument ("positive number, specifying the desired minimum
  distance between nodes. The unit is meter for lonlat data, and in the
  linear unit of the crs for planar data"). Set to `NULL` to not
  densify.

## Value

If `return_poly`, sf, else `out_file`. .parquet mcp written to
`out_file`
