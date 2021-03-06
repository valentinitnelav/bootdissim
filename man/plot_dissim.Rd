% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot_functions.r
\name{plot_dissim}
\alias{plot_dissim}
\title{Plot bootstrapped dissimilarities}
\usage{
plot_dissim(
  metrics_stats,
  size_boot_lines = 0.2,
  alpha_boot_lines = 0.2,
  size_ci = 1,
  linetype_ci = 2,
  size_mean = 1,
  linetype_mean = 1
)
}
\arguments{
\item{metrics_stats}{The output of \code{\link[bootdissim]{get_stats_gg}}.}

\item{size_boot_lines}{Size of the line connecting the bootstrapped values.
Default is 0.2. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.}

\item{alpha_boot_lines}{Alpha parameter of the line connecting the
bootstrapped values. Default is 0.2. Passed to
\code{\link[ggplot2:geom_path]{geom_line}}.}

\item{size_ci}{Size of the line used for the confidence intervals. Default is
\enumerate{
\item Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
}}

\item{linetype_ci}{Type of the line used for the confidence intervals.
Default is 2. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.}

\item{size_mean}{Size of the line connecting the bootstrapped mean values of
the metric. Default is 1. Passed to
\code{\link[ggplot2:geom_path]{geom_line}}.}

\item{linetype_mean}{Type of the line connecting the bootstrapped mean values
of the metric. Default is 1. Passed to
\code{\link[ggplot2:geom_path]{geom_line}}.}
}
\value{
A ggplot object.
}
\description{
Custom plotting function to display the bootstrapped
dissimilarities, the mean line and its confidence intervals.
}
