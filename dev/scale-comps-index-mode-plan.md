# Plan: index-population compositions in `scale_comps()`

Status: implemented 2026-11-23 (steps 1 to 5). Decisions taken: selectivity fixed effects are set to
a reference (default: most frequent level, or the mean for covariates; override with `reference`),
index-mode cells have equal weights unless a `weight` column is given, and roles are given with
explicit `selectivity` and `density` arguments. Averaging over levels of a selectivity term is
not implemented. Recovery tests cover a non-spatial selectivity fixed effect and a spatial fit
with spatially biased sampling. With the TMB backend, either all or none of the random effects can be kept.

## Background

`scale_comps()` currently has two modes, set by `remove_nuisance`.

| Term | Catch weighting (`FALSE`, default) | CPUE mode (`TRUE`) |
|---|---|---|
| `re_form` random effects (e.g. Season, Vessel) | at each stratum's own level | removed (population level) |
| `ffx_form` fixed effects (e.g. Year) | at stratum levels | at stratum levels |
| `cvars` smooths | at the covariate mean in `scale_df` | at the covariate mean in `scale_df` |
| Spatial and spatio-temporal fields | at stratum coordinates | at stratum coordinates |

Catch weighting estimates the size composition of the catch. Each stratum is predicted at its own
levels, including continuous covariates, and the prediction is scaled by the catch of that stratum.
This mode stays as it is.

The CPUE mode is not the right design. It drops random effects only, and the retained-terms split
is arbitrary. The intended meaning is the composition of the abstract population that a CPUE index
tracks, so that the selectivity of the index fleet can be estimated. The prediction should follow
a CPUE analysis.

- Selectivity covariates are dropped (set to a reference).
- In a spatio-temporal analysis, density is predicted over a spatial domain from the spatial and
  spatio-temporal fields and the direct covariates (productivity, SST), after accounting for the
  selectivity covariates.
- In a non-spatial analysis, the index is predicted from Year only, and the length frequencies
  mirror this.
- Predictions are made for each CPUE time step (annual or seasonal), on the same or a similar
  grid and with the same covariates where available.

"Weighting by CPUE" is therefore the wrong description. The output is the composition of the
population described by the index. Some covariates (e.g. SST) can be selectivity or density
terms, and the author has to decide for each term. An earlier version of the code allowed the user
to specify which effects were kept (brms `re_formula`).

## 1. Term roles

Each model term has a role, supplied by the user.

- **Selectivity.** Set to a reference. Examples are vessel, gear, effort-related season, and
  possibly SST.
- **Density.** Kept and predicted. Examples are the spatial and spatio-temporal fields,
  productivity covariates, and possibly SST.
- **Index axis.** The term defining the output time step, usually Year, or Year and Season.

Proposed interface (`fit_model()` is unchanged):

```r
scale_comps(..., mode = "catch")            # current remove_nuisance = FALSE
scale_comps(..., mode = "index",
            selectivity = c("Season", "Vessel", "sst"),   # set to reference
            density      = NULL,                            # default: all other terms
            grid         = pred_grid)                       # spatial models
```

- `remove_nuisance` remains as a deprecated alias (`FALSE` is `mode = "catch"`; `TRUE` is
  `mode = "index"` with `selectivity = fit$re_vars`).
- Terms that are neither listed nor the index axis are not classified automatically. The function
  stops with an error that lists them, so the author makes the decision.

## 2. Reference values for selectivity terms

- Random effects: population level, as now.
- Fixed effects: a reference level (default: most frequent level) or an average over levels. The
  average is closer to a standardised CPUE analysis. **Decision needed.**
- Smooths: a user-supplied value (e.g. the mean), as for `cvars` today.

## 3. Prediction domain

- **Non-spatial:** one prediction per index time step. No catch (`predvar`) is needed. This
  mirrors Year-only index standardisation.
- **Spatial:** predictions on a user-supplied grid that holds coordinates, density covariates and
  time, as in a CPUE spatio-temporal model. Draws are aggregated over cells with weights.
- `grid` replaces `scale_df` and `predvar` in index mode. Catch weighting keeps `scale_df` and
  `predvar`.

## 4. Cell weights

The model is Poisson-multinomial with a `log(n)` offset, so density is not modelled explicitly.
The cell weights in the aggregation therefore need a source.

- (a) Equal weights. The result is the composition of an average cell.
- (b) A `weight` column on the grid, taken from a separate CPUE or abundance model. The length
  frequencies then match the index.

Proposal: implement (b), with (a) as the default. **Decision needed.**

## 5. Output

The long format is unchanged (`.draw`, `bin`, `tot_by_bin`, plus the index axis), so
`scaled_comp_plot()` and `scaled_ridge_plot()` work as they are. The selectivity of the index fleet
is estimated in the stock assessment model, which takes this composition as input. No helper
for it is needed in this package.

## 6. Work breakdown

1. Add `mode`, `selectivity`, `density` and `grid`, with validation and the deprecation alias.
2. Rewrite the sdmTMB and brms prediction paths. brms uses `re_formula` plus reference values.
   sdmTMB uses reference values for fixed terms and `re_form_iid = NA` for random effects.
3. Add grid aggregation with weights.
4. Add a recovery test. The simulation has a known population composition and sampling biased by
   vessel and location.
5. Update `?scale_comps`, `NEWS.md` and the vignettes. The SST/ENSO section of
   `spatial_bias_simulation.Rmd` becomes the index-mode example.

## Open decisions

1. Fixed-effect selectivity terms: reference level (most frequent) or average over levels?
2. Weights in index mode: equal cells by default, plus an optional `weight` column from a separate
   CPUE model?
3. Interface: explicit `selectivity`/`density` arguments with an error on unclassified terms, or
   roles inferred from `re_form` membership?
