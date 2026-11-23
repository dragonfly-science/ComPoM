# Plot the step results

One panel per step, with the statistic by level of `by` for the current
step (orange, with its 95% interval), the previous step (skyblue) and
earlier steps (grey). Rows are the components (catch and index).

## Usage

``` r
plot_step_comps(step, xlab = step$by, ylab = step$stat_label)
```

## Arguments

- step:

  Output of
  [`step_comps()`](https://dragonfly-science.github.io/ComPoM/reference/step_comps.md).

- xlab, ylab:

  Axis labels. `ylab` defaults to the statistic.

## Value

A ggplot.

## See also

[`step_comps()`](https://dragonfly-science.github.io/ComPoM/reference/step_comps.md)
