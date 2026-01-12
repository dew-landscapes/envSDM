# Lookup for metrics relevant to evaluating SDMs

Lookup for metrics relevant to evaluating SDMs

## Usage

``` r
sdm_metrics
```

## Format

A data frame with 18 rows and 5 variables:

- metric:

  Character. Short name of metric used in outputs

- high_good:

  Logical. Do high values of metric correlate with a good model?

- is_thresh:

  Logical. Is this metric dependent on setting a threshold?

- summary_mets:

  Logical. Which metrics to use in summarising?

- within_mets:

  Logical. Included here for compatibility with
  [`envFunc::make_metric_df()`](https://rdrr.io/pkg/envFunc/man/make_metric_df.html)

- across_mets:

  Logical. Included here for compatibility with
  [`envFunc::make_metric_df()`](https://rdrr.io/pkg/envFunc/man/make_metric_df.html)
