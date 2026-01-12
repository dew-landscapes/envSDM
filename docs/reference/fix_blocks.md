# Make sure all blocks achieve min_fold_n presences

Make sure all blocks achieve min_fold_n presences

## Usage

``` r
fix_blocks(blocks, pres, min_fold_n = 8, pres_val = 1)
```

## Arguments

- blocks:

  Vector of block ids

- pres:

  Vector of presence/absence data

- min_fold_n:

  Minimum number of presences in a block

- pres_val:

  Value in `pres` that represent presence

## Value

Vector of adjusted block ids ensuring each block achieves min_fold_n
presences
