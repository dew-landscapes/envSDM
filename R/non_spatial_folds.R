#' Generate non spatial folds
#'
#' @param use_folds Numeric. Total number of folds
#' @param data Dataframe with `pa_col`
#' @param pa_col Name of column in data containing p/a values
#' @param pres_val Value in `pa_col` identifying presences
#' @param min_in_fold `min_fold_n`
#' @param max_attempts How many attempts to make to achieve `min_in_fold`
#' presences within each fold?
#'
#' @returns
#' @export
#' @keywords internal
#'
#' @examples
non_spatial_folds <- function(use_folds
                              , data
                              , pa_col = "pa"
                              , pres_val = 1
                              , min_in_fold = 5
                              , max_attempts = 99
                              ) {

  flds <- function(fs = 1:use_folds, val = 1) {

    n_folds <- length(unique(fs))

    sample(fs
           , sum(data[,pa_col] == val)
           , replace = TRUE
           , prob = rep(1 / n_folds, n_folds)
           )

  }

  p_fold <- flds()
  counter <- 0

  # attempt to get min_in_fold presences in each fold
  while(any(min(table(p_fold)) < min_in_fold, counter <= max_attempts)) {

    counter <- counter + 1
    p_fold <- flds()

  }

  # if that fails, use fix_folds
  if(min(table(p_fold)) < min_in_fold) {

    p_fold <- fix_folds(folds = p_fold
                        , pres = rep(1, length(p_fold))
                        , min_fold_n = min_in_fold
                        , pres_val = pres_val
                        )

  }

  c(p_fold
    , flds(fs = unique(p_fold)
           , val = 0
           )
    )

}
