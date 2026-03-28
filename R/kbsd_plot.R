#' Plotting the KBSD
#'
#' @param result_df resulting data.frame of kbsd()
#'
#' @return plot
#' @importFrom stats sd quantile
#' @importFrom graphics abline bxp points
#' @export
kbsd_plot <- function(result_df){
  groups <- split(result_df$diagnostic, result_df$intervention)

  # Quantile berechnen
  q <- lapply(groups, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)

  # In Matrix umwandeln
  stats <- sapply(q, identity)

  # Boxplot zeichnen
  bx <- bxp(list(
    stats = stats,
    n = sapply(groups, length),
    names = names(groups)
  ),
  ylim = c(0, max(result_df$diagnostic, na.rm = TRUE)),
  xlab = "intervention",
  ylab = "EDP",
  main = "KBSD: Adjusted Boxplots"
  )
  abline(h = 0, col = "lightgray", lwd = 1, lty = 2)

  # Ausreißer (unter 5% oder über 95%) hinzufügen
  for(i in seq_along(groups)) {
    x <- groups[[i]]
    out <- x[x < stats[1,i] | x > stats[5,i]]
    points(rep(i, length(out)), out, pch = 16, cex=0.5)
  }
}
