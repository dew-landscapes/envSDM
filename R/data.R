
#' Lookup for metrics relevant to evaluating SDMs
#'
#'
#' @format A data frame with `r nrow(sdm_metrics)` rows and
#' `r ncol(sdm_metrics)` variables:
#' \describe{
#'   \item{metric}{Character. Short name of metric used in outputs}
#'   \item{high_good}{Logical. Do high values of metric correlate with a good model?}
#'   \item{is_thresh}{Logical. Is this metric dependent on setting a threshold?}
#'   \item{summary_mets}{Logical. Which metrics to use in summarising?}
#'   \item{within_mets}{Logical. Included here for compatibility with
#'   `envFunc::make_metric_df()`}
#'   \item{across_mets}{Logical. Included here for compatibility with
#'   `envFunc::make_metric_df()`}
#'   ...
#' }
"sdm_metrics"
