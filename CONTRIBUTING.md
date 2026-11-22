# Contributing to ComPoM

Bug reports and pull requests are welcome at <https://github.com/dragonfly-science/ComPoM>.

## Bug reports

Open an issue with a minimal reproducible example. `ex_comp` and `ex_catch` (see `?ex_data`)
are bundled so that examples can use them. Include `sessionInfo()`, and the sdmTMB and brms
versions.

## Development

```r
devtools::load_all()
devtools::document()                      # regenerate man/ and NAMESPACE (roxygen2)
devtools::test()                          # unit and regression tests
devtools::check(vignettes = FALSE)        # as in CI
```

* Edit the roxygen comments in `R/`, not the files in `man/`. CI fails if `man/` or `NAMESPACE`
  are out of date.
* Add a test for each change in behaviour, and a line to `NEWS.md`.
* The brms tests need `rstan` and a C++ toolchain, and are skipped if `rstan` is not installed.
* The vignettes are slow (the spatial vignette takes about 10 minutes). Render them with
  `rmarkdown::render()` from the `vignettes/` directory when changing the model or the plots.
* `data-raw/ex_data.R` regenerates the example data.
