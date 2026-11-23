# Coefficient, distribution and influence plot for one term

A coefficient-distribution-influence plot (Bentley et al. 2011) for one
term of
[`influ_comps()`](https://dragonfly-science.github.io/ComPoM/reference/influ_comps.md),
with the axes shared between the panels.

- Centre: the distribution of the data. Bubbles show the share of the
  data in each level of the term (x) in each year (y).

- Top, sharing the level axis (x): the effect of each level on every bin
  (y), as the log of the multiplier relative to the geometric mean over
  levels. Bars extend to the right (more fish in the bin than average)
  or to the left (fewer), with 95% intervals. Bar lengths are scaled to
  the largest effect, and a scale bar gives the reference length. With
  many levels, a heat map replaces the bars. Cells whose interval
  includes zero are faded.

- Right, sharing the year axis (y): the influence of the term on the
  annual composition in every bin (x), as the log ratio of the predicted
  proportion with and without the term. Cells whose interval includes
  zero are faded.

Orange is below and skyblue above the reference.

## Usage

``` r
plot_cdi_comps(
  influ,
  term,
  component = NULL,
  var = NULL,
  style = c("auto", "bars", "tile"),
  max_levels = 10,
  xlab = "Length bin (cm)"
)
```

## Arguments

- influ:

  Output of
  [`influ_comps()`](https://dragonfly-science.github.io/ComPoM/reference/influ_comps.md).

- term:

  Term (a name in `influ$groups`, i.e. a variable, or the joint
  random-effect label).

- component:

  `"catch"` or `"index"`, for the distribution and influence panels.

- var:

  Variable whose levels are shown, for a term that groups several
  variables. Default: the first.

- style:

  `"auto"` (bars for up to `max_levels` levels, otherwise a heat map),
  `"bars"` or `"tile"`.

- max_levels:

  Maximum number of levels shown as bars with `style = "auto"`.

- xlab:

  Label of the bin axis.

## Value

A patchwork of ggplots.

## References

Bentley, N., Kendrick, T.H., Starr, P.J., Breen, P.A. (2012). Influence
plots and metrics: tools for better understanding fisheries
catch-per-unit-effort standardizations. ICES Journal of Marine Science
69, 84-88.

## See also

[`influ_comps()`](https://dragonfly-science.github.io/ComPoM/reference/influ_comps.md)
