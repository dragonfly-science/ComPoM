% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_proc.R
\name{scaled_ridge_plot}
\alias{scaled_ridge_plot}
\title{Plot scaled compositions as ridgeline plot}
\usage{
scaled_ridge_plot(
  scaled_comp = NULL,
  grps = NULL,
  scales = "free_y",
  plot_y = NULL,
  plot_facet = NULL
)
}
\arguments{
\item{scaled_comp}{data frame used for scaling}

\item{grps}{groups to do scaling}

\item{scales}{free y axis?}

\item{plot_y}{y axis for ridgeline}

\item{plot_facet}{faceting ridgelines (e.g., monthly or quarterly ridges within years)}
}
\description{
Plot scaled compositions as ridgeline plot
}
