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

  flds <- function(fs = use_folds, val = 1) {

    sample(x = fs
           , size = sum(data[,pa_col] == val)
           , replace = TRUE
           , prob = rep(1 / fs, fs)
           )

  }

  result <- flds()
  counter <- 0

  # attempt to get min_in_fold presences in each fold
  while(all(min(table(result)) < min_in_fold, counter <= max_attempts)) {

    counter <- counter + 1
    result <- flds()

  }

  # if that fails, use fix_folds
  if(min(table(result)) < min_in_fold) {

    result <- fix_folds(folds = result
                        , pres = rep(1, length(result))
                        , min_fold_n = min_in_fold
                        , pres_val = pres_val
                        )

  }

  if(sum(data[[pa_col]] != pres_val)) {

    result <- c(result
                , flds(fs = use_folds
                       , val = 0
                       )
                )

  }

  return(result)

}
