# Statistics of a composition

Helpers for the `stat` argument of
[`step_comps()`](https://dragonfly-science.github.io/ComPoM/reference/step_comps.md)
and
[`influ_comps()`](https://dragonfly-science.github.io/ComPoM/reference/influ_comps.md).
A statistic is a function of the bin values and the proportions in each
bin that returns one number.

## Usage

``` r
comp_quantile(prob)
```

## Arguments

- prob:

  Probability of the quantile, between 0 and 1.

## Value

A function with arguments `bin` and `p`.

## Examples

``` r
comp_quantile(0.9)(bin = c(10, 20, 30), p = c(0.2, 0.5, 0.3))
#> [1] 30
```
