# Parse a formula

Constructs the Poisson-multinomial model formula, separating fixed
effects (variables of primary interest, entered as `factor(var):bin`
interactions) from nuisance random effects (entered as `(1|bin:var)` IID
random intercepts). Called internally by
[`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md);
it is exported to allow inspection of the formula and model data before
fitting.

## Usage

``` r
parse_form(
  data,
  backend = "brms",
  ffx_form = "Year",
  re_form = NULL,
  cvars = NULL,
  spl = "2D",
  knots = "c(5,10)",
  coords = NULL,
  time = NULL
)
```

## Arguments

- data:

  Prepared data frame from
  [`data_prep()`](https://dragonfly-science.github.io/ComPoM/reference/data_prep.md)

- backend:

  the backend used for model fitting (`"brms"` or `"TMB"`)

- ffx_form:

  Character string of primary fixed-effect variables separated by `+`,
  e.g. `"Year"`. These become `factor(Year):bin` fixed-effect
  interactions in the model.

- re_form:

  Character string of nuisance random-effect variables separated by `+`,
  e.g. `"Season + Gear"`. These become `(1|bin:Season)` IID random
  intercepts.
  [`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md)
  evaluates them at each stratum's level (`mode = "catch"`) or sets them
  to the population level (`mode = "index"`).

- cvars:

  Character vector of continuous covariates to be fitted as smooths with
  bin (see Details). `NULL` = none.

- spl:

  Smooth type: `"2D"`, `"by_factor"` or `"by"` (see Details).

- knots:

  Basis dimension of the smooths, as a character string: a single number
  (e.g. `"4"`) for `spl = "by_factor"` or `"by"`, or `"c(5,10)"`
  (covariate, bin) for `spl = "2D"`.

- coords:

  Character vector of coordinate column names that are kept numeric
  (spatial TMB models only)

- time:

  Name of the time column that is kept numeric/unchanged (TMB backend
  only)

## Value

A list with the model formula (`form`), the model data with grouping
variables converted to factors (`data`), and the names of the nuisance
random-effect variables (`re_vars`).

## Details

The linear predictor always contains bin-specific intercepts
(`0 + bin`). Continuous covariates in `cvars` are added as smooths
according to `spl`:

- `"2D"`:

  `t2(cvar, cbin, k = knots)`, a tensor-product smooth over the
  covariate and the numeric bin value. Requires numeric bins (e.g.
  length). With the `"TMB"` backend this smooth gave a
  non-positive-definite Hessian in all simulations tried, because its
  marginal bin component is collinear with the bin-specific intercepts
  (the numeric bin value has only as many distinct values as there are
  bins). Use `"by_factor"` with that backend.

- `"by_factor"`:

  `s(cvar, by = bin, k = knots)`, a separate 1D smooth of the covariate
  for each bin. Also applicable to non-ordered categories (e.g.
  species).

- `"by"`:

  Legacy option. Separate `t2()` smooths for each non-reference bin via
  bin dummy columns.

## Examples

``` r
dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
pf <- parse_form(dat, backend = "TMB", ffx_form = "Year", re_form = "Season + Vessel")
#> Model formula: tot_by_bin ~ 0 + bin + factor(Year):bin + (1 | bin:Season) + (1 | bin:Vessel)
pf$form
#> tot_by_bin ~ 0 + bin + factor(Year):bin + (1 | bin:Season) + 
#>     (1 | bin:Vessel)
#> <environment: 0x5fea16d1b2b0>
pf$re_vars
#> [1] "Season" "Vessel"
```
