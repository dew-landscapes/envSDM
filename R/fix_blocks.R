#' Make sure all blocks achieve min_fold_n presences
#'
#' @param blocks Vector of block ids
#' @param pres Vector of presence/absence data
#' @param min_fold_n Minimum number of presences in a block
#' @param pres_val Value in `pres` that represent presence
#'
#' @returns Vector of adjusted block ids ensuring each block achieves min_fold_n
#' presences
#' @export
#' @keywords internal
#' @examples
fix_blocks <- function(blocks, pres, min_fold_n = 8, pres_val = 1) {

  k_folds <- length(unique(blocks))

  blocks_p <- blocks[pres == pres_val]

  if(any(table(blocks_p) < min_fold_n, length(setdiff(1:k_folds, unique(blocks_p))) > 0)) {

    how_many_below_thresh <- sum(purrr::map_dbl(1:k_folds, \(x) sum(blocks_p == x)) < min_fold_n)

    old_k_folds <- k_folds

    k_folds <- old_k_folds - how_many_below_thresh

    blocks_adj <- tibble::tibble(fold_ids = 1:old_k_folds) %>%
      dplyr::mutate(n = purrr::map_dbl(fold_ids, \(x) sum(blocks_p == x))) %>%
      dplyr::mutate(fold_ids_adj = forcats::fct_lump_n(as.factor(fold_ids), k_folds - 1, w = n)) %>%
      dplyr::mutate(fold_ids_adj = as.numeric(fold_ids_adj)) %>%
      dplyr::distinct()

    blocks <- blocks_adj$fold_ids_adj[match(blocks, blocks_adj$fold_ids)]

  }

  return(blocks)
}
