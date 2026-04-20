# Generate non spatial folds

Generate non spatial folds

## Usage

``` r
non_spatial_folds(
  use_folds,
  data,
  pa_col = "pa",
  pres_val = 1,
  min_in_fold = 5,
  max_attempts = 99
)
```

## Arguments

- use_folds:

  Numeric. Total number of folds

- data:

  Dataframe with `pa_col`

- pa_col:

  Name of column in data containing p/a values

- pres_val:

  Value in `pa_col` identifying presences

- min_in_fold:

  `min_fold_n`

- max_attempts:

  How many attempts to make to achieve `min_in_fold` presences within
  each fold?
