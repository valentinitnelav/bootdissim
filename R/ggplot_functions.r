#' @title Prepare bootstrapped results for `ggplot` friendly data format.
#'
#' @description Helper function to convert the bootstrapped dissimilarities from
#'   matrix format to `ggplot2` friendly data format.
#'
#' @param x A matrix of bootstrapped dissimilarities from
#'   \code{\link[bootdissim]{boot_dissim}}
#'
#' @param probs A numeric vector of two probabilities in the interval `[0, 1]`.
#'   Used for building the lower and upper bounds of the confidence intervals.
#'   Defaults to `c(0.025, 0.975)`, which corresponds to a 95\\% confidence
#'   interval.
#'
#' @return A list of two data frames to be used with
#'   \code{\link[ggplot2]{ggplot}}. First one, `stats_df` is a 4 columns data
#'   frame, containing the sample size/number of sampled interactions
#'   (`spl_size`) at each step and the corresponding mean dissimilarity (`mean`)
#'   together with its bootstrap confidence interval limits (`ci_low`, `ci_up`).
#'   The second one, `lines_df`, contains all the bootstrapped dissimilarities
#'   (`value`) at each bootstrap iteration (`simulation_id`) and its
#'   corresponding sample size/sampled interactions (`spl_size`). Can be used
#'   for enhancing visual effect when plotting the mean bootstrap values.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table melt setDT
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom matrixStats rowMeans2 rowQuantiles
#'
#' @export
#'
#' @md
get_stats_gg <- function(x, probs = c(0.025, 0.975)){
  if (! is.matrix(x) ) stop("Expecting a matrix")

  stats_df <- data.frame(spl_size = as.integer(rownames(x))) %>%
    dplyr::mutate(mean   = matrixStats::rowMeans2(x, na.rm = TRUE),
                  ci_low = matrixStats::rowQuantiles(x, probs = probs[1], na.rm = TRUE) %>% unname,
                  ci_up  = matrixStats::rowQuantiles(x, probs = probs[2], na.rm = TRUE) %>% unname)

  lines_df <- x %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "spl_size") %>%
    setDT() %>%
    data.table::melt(id.vars = "spl_size",
                     variable.name = "simulation_id") %>%
    dplyr::mutate(spl_size = as.integer(spl_size),
                  simulation_id = as.integer(simulation_id))

  return(list(stats_df = stats_df,
              lines_df = lines_df))
}


#' @title Plot bootstrapped dissimilarities
#'
#' @description Custom plotting function to display the bootstrapped
#'   dissimilarities, the mean line and its confidence intervals.
#'
#' @param metrics_stats The output of \code{\link[bootdissim]{get_stats_gg}}.
#'
#' @param size_boot_lines Size of the line connecting the bootstrapped values.
#'   Default is 0.2. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param alpha_boot_lines Alpha parameter of the line connecting the
#'   bootstrapped values. Default is 0.2. Passed to
#'   \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param size_ci Size of the line used for the confidence intervals. Default is
#'   1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param linetype_ci Type of the line used for the confidence intervals.
#'   Default is 2. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param size_mean Size of the line connecting the bootstrapped mean values of
#'   the metric. Default is 1. Passed to
#'   \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param linetype_mean Type of the line connecting the bootstrapped mean values
#'   of the metric. Default is 1. Passed to
#'   \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#'
#' @export
#'
#' @md
plot_dissim <- function(metrics_stats,
                        size_boot_lines = 0.2,
                        alpha_boot_lines = 0.2,
                        size_ci = 1,
                        linetype_ci = 2,
                        size_mean = 1,
                        linetype_mean = 1){

  stats_df <- metrics_stats[["stats_df"]]
  lines_df <- metrics_stats[["lines_df"]]

  plots_gg <- ggplot() +
    # Bootstrap lines
    geom_line(data = lines_df,
              aes(x = spl_size,
                  y = value,
                  group = simulation_id),
              size = size_boot_lines,
              alpha = alpha_boot_lines) +
    # Lower confidence interval bound
    geom_line(data = stats_df,
              aes(x = spl_size,
                  y = ci_low),
              size = size_ci,
              linetype = linetype_ci) +
    # Upper confidence interval bound
    geom_line(data = stats_df,
              aes(x = spl_size,
                  y = ci_up),
              size = size_ci,
              linetype = linetype_ci) +
    # Average line
    geom_line(data = stats_df,
              aes(x = spl_size,
                  y = mean),
              size = size_mean,
              linetype = linetype_mean)

  return(plots_gg)
}
