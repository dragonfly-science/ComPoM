# Geometric mean; used by Fx_plot() to express bin proportions relative to
# their geometric mean across levels of the plotted grouping variable.
gmean <- function(x) exp(mean(log(x)))
