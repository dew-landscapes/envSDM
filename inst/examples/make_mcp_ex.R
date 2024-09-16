
# prep -------
use_crs = 8059
latlong_crs = 4326

use_clip <- sa_vect |>
  terra::unwrap() |>
  sf::st_as_sf() |>
  dplyr::summarise() |>
  sf::st_transform(crs = use_crs) |>
  sf::st_make_valid()

out_dir <- file.path(system.file(package = "envSDM"), "examples")

# mcp -------

mcp_prep <- clean_end |>
  dplyr::distinct(long, lat, taxa) |>
  tidyr::nest(data = -c(taxa)) |>
  dplyr::filter(purrr::map_lgl(data, \(x) nrow(x) > 2)) |>
  dplyr::mutate(out_file = fs::path(out_dir, taxa, "mcp.parquet"))


purrr::walk2(mcp_prep$data
             , mcp_prep$out_file
             , \(x, y) make_mcp(presence = x
                                , out_file = y
                                , force_new = FALSE
                                , in_crs = latlong_crs
                                , out_crs = use_crs
                                , buf = 0
                                , clip = use_clip
             )
)
