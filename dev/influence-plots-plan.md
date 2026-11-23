# Plan: step and CDI plots for compositions

Status: proposal, nothing implemented. Written 2026-11-23.

## What exists

Source: `~/dragonfly/Kahawai/generalised_influ` (Ginflu 0.1.1) and
`~/dragonfly/Inshore_CPUE/influ2`. The comparison of the two is in `~/dragonfly/Kahawai/influ2_issue.md`.

- **Step plot** (`get_step()`, `plot_step()`). The model is refitted with terms added one at a
  time (intercept, then each term, then the spatial and spatio-temporal fields for sdmTMB). The
  annual index of each refit is stored (`Model`, `Index`, `level`, `mean`, `Lower`, `Upper`) with
  deviance and AIC changes per step. The plot shows one panel per step, with the current,
  previous and earlier indices overlaid.
- **CDI plot** (`cdi_plot_with_indicators()`). One figure per term with three panes. The
  coefficient pane shows the effect of each level of the term. The distribution pane shows the share
  of the data in each level by year. The influence pane shows the contribution of the term to the
  annual index (Bentley et al. 2011), with summary metrics.
- Both work on a scalar response per observation, so one coefficient per level and one index per
  year. A ComPoM model has one set of coefficients per bin (`factor(var):bin`, `(1|bin:var)`), and
  its output is a composition, not a scalar.

## Approach

Reduce the composition to a scalar summary S per year, and compute everything with the existing
prediction machinery (`scale_comps()`), so that the same code handles brms and sdmTMB fits and all
term types (fixed effects, random effects, covariate smooths, spatial fields).

- S is the mean bin (mean length) by default. Options are a proportion above or below a bin, or any
  user function of the composition (`stat = function(bin, p) ...`).
- The full composition per year can be shown as small multiples (one panel per year, bins on the
  x axis) in place of S.

### Step plot

1. Refit with terms added in the order given (default: `ffx_form` terms, then `re_form` terms,
   then covariate smooths, then spatial fields). The first step is the bin-only model (one overall
   composition).
2. For each refit, predict the annual composition and compute S per year with uncertainty from the
   draws.
3. Output: a long data frame in the Ginflu step format (`Model`, `level`, `mean`, `Lower`,
   `Upper`), plus a summary table (df, logLik, AIC, deviance change where available), so that
   `Ginflu::plot_step()` can plot it. A ComPoM plot function with `theme_compom()` is the
   alternative.
4. The two components are two settings of the prediction at each step.
   - **Catch component:** `mode = "catch"`. Each stratum is predicted at its own levels with the
     terms present and weighted by catch.
   - **Index component:** `mode = "index"` on the grid, with the terms present. Selectivity terms
     are set to a reference once they have been added, so the step shows when they stop affecting
     the index.

Cost: one refit per term. TMB fits take seconds. brms fits take minutes each, so brms steps need an
option to reuse fewer chains or iterations.

### CDI plot

For a term k and year y, the influence is the difference in S between the prediction with the term
at the levels observed in year y and the prediction with the term at its reference level (or mean,
for a smooth). This is the same quantity as the Bentley influence (the average contribution of the
term over the observations of the year), computed on the composition scale.

- **Coefficient pane (bin-resolved).** The effect of each level of the term on every bin, not on
  S. The effect of level l on bin b is the log ratio of the predicted proportion at level l to
  the proportion at the reference level (or to the geometric mean across levels), with intervals
  from the draws. This is the multiplier that `Fx_plot()` already shows, reorganised.
  - Bins run vertically (the length axis). For each level, a bar starts at zero and extends to the
    left (fewer fish in the bin than at the reference) or to the right (more), filled with a
    diverging scale (orange below zero, skyblue above).
  - Few levels (about 8 or fewer, e.g. season, year): one panel per level, side by side.
  - Many levels (e.g. 40 vessels): a heat map with levels on x, bins on y and the effect as a
    diverging fill, ordered by the mean effect.
  - Proportions sum to one, so the effects are relative. One bin going up means others go down.
    The log ratio is therefore centred across bins (centred log ratio) when the reference is the
    geometric mean, so that the bars sum to zero.
  - Covariate smooths: the same bars at a few values of the covariate (quartiles), or a heat map of
    covariate value by bin. Spatial fields: one map per bin, as in `spatial_field_plot()`.
- **Distribution pane.** The share of the samples (or of the catch in catch mode) in each level by
  year. This needs no model, only the data and the strata.
- **Influence pane.** The influence by year, as above, with the same summary metrics as Ginflu.
  The influence is also bin-resolved: a heat map with year on x, bin on y and the log ratio of
  the predicted proportion with and without the term as a diverging fill. S is shown as an
  optional line beneath it, because Ginflu's metrics (overall influence, trend) need a scalar.
- **Catch against index.** The influence of a term on the catch composition, and on the index
  composition, can be shown side by side. For selectivity terms, the influence on the index is
  zero by construction. The difference between the two is the effect of the term on the sampling,
  which is the quantity the stock assessment needs to interpret the index fleet.
- Continuous covariates and spatial fields are binned the way Ginflu bins numeric terms
  (quantile classes for the first term, pretty breaks for the second). For spatial fields, the
  coefficient pane becomes a map of the field's effect on S, and the distribution pane a map of the
  sample density.

## Where to put the code

1. **Compute in ComPoM**, return data frames in the Ginflu long format, and plot with Ginflu or
   with ComPoM functions that use `theme_compom()`. This keeps one definition of the predictions.
2. **Contribute to Ginflu/influ2.** The calculations are specific to compositions, so they only
   suit the packages if they accept a user-supplied prediction function. That fits the plan in
   `influ2_issue.md` to put the calculations in influ2.

Recommendation: option 1 first, with the data frame format kept compatible, then decide about
option 2 after the results look right.

## Work breakdown

1. `comp_summary()`: S per year and draw from a `scale_comps()` result, with `stat` options.
2. `step_comps()`: refits and S per step for both modes. Tests on the simulated vessel data, where
   the true effect of each term is known.
3. `influ_comps()`: influence and coefficient tables per term and year. Tests: the influence of a
   term with no effect is zero, and the influences sum approximately to the difference between the
   full prediction and the bin-only prediction.
4. `plot_step_comps()` and `plot_cdi_comps()` with `theme_compom()`, or adapters to the Ginflu
   plot functions.
5. A vignette section on the vessel simulation.

## Decisions needed

1. Summary S for the step plot and the influence metrics: mean length, or something else (for
   example proportion above a threshold). The coefficient and influence panes are bin-resolved and
   do not need S.
2. Plotting: new ComPoM functions, or adapters to Ginflu so that the plots match the existing
   CPUE reports?
3. Whether the step order should follow the model formula or be set by the user.
4. Whether brms models are in scope for the first version, given the refit cost.
