#' creates a half distance vector according to rule-of-thumb
#'
#' @param data Observed data.frame
#' @param interventions a vector of strings naming all intervention variables
#' @param remove a vector of strings naming dimensions not to be considered, e.g. the outcome variable
#' @param no_extrapolation a vector of strings naming all dimensions where extrapolation is excluded, e.g. binary variables
#'
#' @return a named vector
#' @export
get_disthalf_vec <- function(data, interventions = NULL, remove = NULL, no_extrapolation = NULL){

  # standard deviations (numeric only)
  res <- sapply(data, function(x) {
    if (is.factor(x) || is.logical(x)) {
      0
    } else {
      sd(x, na.rm = TRUE)
    }
  })

  # no extrapolation, e.g. binary variables
  if (!is.null(no_extrapolation)) {
    res[names(res) %in% no_extrapolation] <- 0
  }

  # remove outcome etc.
  if (!is.null(remove)) {
    res <- res[setdiff(names(res), remove)]
  }

  # halve intervention variables
  if (!is.null(interventions)) {
    res[interventions] <- res[interventions] / 2
  }

  return(res)
}
