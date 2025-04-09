#' Generate non spatial blocks
#'
#' @param use_folds Numeric. Total number of folds
#' @param data Dataframe with `pa_col`
#' @param pa_col Name of column in data containing p/a values
#' @param pres_val Value in `pa_col` identifying presences
#'
#' @returns
#' @export
#' @keywords internal
#'
#' @examples
non_spatial_blocks <- function(use_folds, data, pa_col = "pa", pres_val = 1) {

  c(sample(1:use_folds
           , sum(data[,pa_col] == 1)
           , replace = TRUE
           , prob = rep(1 / use_folds, use_folds)
           )
    , sample(1:use_folds
             , sum(data[,pa_col] == 0)
             , replace = TRUE
             , prob = rep(1 / use_folds, use_folds)
             )
    )

}
