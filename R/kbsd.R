library(Rfast)


# Gaussian (rbf) kernel
.gaussian_kernel <- function(x, point, disthalf) {
  gamma <- -log(0.5)/(disthalf)^2
    exp(-gamma * (x-point)^2 )
}

# single observation diagnostic
.compute_single_diagnostic <- function(data_observed, obs_intervened, disthalf_vec, kernel = .gaussian_kernel){
  cols_used <- intersect(colnames(data_observed), names(disthalf_vec))
  # pre-extract bandwidths in correct order
  disthalf <- disthalf_vec[cols_used]
  # vectorized row-wise computation
  res_total <- sum(apply(data_observed, 1, function(row_i){
    if (length(cols_used) == 0) return(1)
    prod(
      mapply(
        function(colname, h) {
          kernel(
            row_i[[colname]],
            obs_intervened[[colname]],
            disthalf = unname(h)
          )
        },
        cols_used,
        disthalf
      )
    )
  }))
  res_total
}

# ----------------------------
# Main diagnostic function
# ----------------------------

#' Kernel-based sparsity diagnostic
#'
#' @param data Observed data.frame
#' @param int_data_list List of intervened on datasets
#' @param disthalf_vec Named numeric vector of half-distances
#' @param type Either "Rfast" (fast) or "base"
#' @param kernel Kernel function or "gaussian"
#' @param parallel Logical
#'
#' @return data.frame with diagnostics
#' @export
kbsd <- function(data, int_data_list, disthalf_vec, type="Rfast", kernel="gaussian", parallel=TRUE){

  # to handle division by zero issues
  disthalf_vec[disthalf_vec == 0] <- .Machine$double.xmin

  # set the kernel
  if (kernel=="gaussian"){
    kernel_used <- .gaussian_kernel
  } else {
    kernel_used <- kernel
  }

  # only keep relevant variables for the observed data
  original_data <- data
  vartoconsider <- names(disthalf_vec)
  data <- data[vartoconsider]

  # apply the diagnostic to all interventions
  result_df <- data.frame(matrix(NA, nrow=0, ncol=3))
  names(result_df) <- c("intervention", "diagnostic", "observation")
  for (i in seq_along(int_data_list)){
    # intervened on data
    data_intervened <- int_data_list[[i]][vartoconsider]

    # calculation of the diagnostic
    if (type=="Rfast") {
      # make sure correct Rfast version is used
      if (utils::packageVersion("Rfast") > "2.1.5.1") {
        warning(
          'Rfast > 2.1.5.1 is known to be buggy and leads to incorrect results for this function.
To install another version, for example 2.1.5.1, restart R (and RStudio if applicable) and run:
remotes::install_version("Rfast", version = "2.1.5.1")
Check which version is in use with:
utils::packageVersion("Rfast")',
          call. = FALSE
        )
      }
      # use Rfast to calculate the diagnostic
      dat_int_scaled <- sweep(data_intervened, 2, 1/disthalf_vec, FUN='*')
      dat_scaled <- sweep(data, 2, 1/disthalf_vec, FUN='*')
      dia <- Rfast::rowsums(exp(log(0.5)*Rfast::dista(dat_int_scaled, dat_scaled, square = TRUE)), parallel = parallel)
    } else {
      dia <- apply(data_intervened, 1, function(x){.compute_single_diagnostic(data, x, disthalf_vec,
                                                                    kernel=kernel_used)})
    }
    new_data <- data.frame(intervention= rep(i, length(dia)), diagnostic=unname(dia), observation=1:length(dia) )
    result_df <- rbind(result_df, new_data)
  }

  # return result
  return(result_df)
}
