# ComPoM plot theme

The `cowplot` theme with faded major gridlines, used by all plotting
functions in the package. Use it (or
`ggplot2::theme_set(theme_compom())`) to match other figures to the
package plots.

## Usage

``` r
theme_compom(font_size = 11)
```

## Arguments

- font_size:

  Base font size.

## Value

A ggplot2 theme.

## Examples

``` r
ggplot2::ggplot(ex_comp, ggplot2::aes(length, count)) +
  ggplot2::geom_col() +
  theme_compom()
```
