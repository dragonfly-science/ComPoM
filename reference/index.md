# Package index

## Workflow

Prepare data, fit the model and scale the predictions to catch.

- [`data_prep()`](https://dragonfly-science.github.io/ComPoM/reference/data_prep.md)
  : Prepare compositional count data
- [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)
  : Fit the Poisson-multinomial model
- [`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md)
  : Scale predicted compositions to catch or to an index population
- [`parse_form()`](https://dragonfly-science.github.io/ComPoM/reference/parse_form.md)
  : Parse a formula

## Diagnostics and effects

- [`post_pred_group()`](https://dragonfly-science.github.io/ComPoM/reference/post_pred_group.md)
  : Posterior predictive check by group
- [`Fx_plot()`](https://dragonfly-science.github.io/ComPoM/reference/Fx_plot.md)
  : Plot bin-specific effects as multipliers
- [`cvar_smooth_plot()`](https://dragonfly-science.github.io/ComPoM/reference/cvar_smooth_plot.md)
  : Plot marginal smooth effects of a continuous covariate by bin
- [`spatial_field_plot()`](https://dragonfly-science.github.io/ComPoM/reference/spatial_field_plot.md)
  : Plot predicted spatial field by bin

## Influence of model terms

Step and coefficient-distribution-influence plots for compositions.

- [`step_comps()`](https://dragonfly-science.github.io/ComPoM/reference/step_comps.md)
  : Add model terms one at a time and track a statistic of the
  composition
- [`plot_step_comps()`](https://dragonfly-science.github.io/ComPoM/reference/plot_step_comps.md)
  : Plot the step results
- [`influ_comps()`](https://dragonfly-science.github.io/ComPoM/reference/influ_comps.md)
  : Influence of model terms on the composition
- [`plot_cdi_comps()`](https://dragonfly-science.github.io/ComPoM/reference/plot_cdi_comps.md)
  : Coefficient, distribution and influence plot for one term
- [`comp_quantile()`](https://dragonfly-science.github.io/ComPoM/reference/comp_quantile.md)
  : Statistics of a composition

## Plot theme

- [`theme_compom()`](https://dragonfly-science.github.io/ComPoM/reference/theme_compom.md)
  : ComPoM plot theme

## Scaled compositions

- [`scaled_comp_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_comp_plot.md)
  : Plot scaled compositions
- [`scaled_ridge_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_ridge_plot.md)
  : Plot scaled compositions as ridgeline plot

## Data

- [`ex_data`](https://dragonfly-science.github.io/ComPoM/reference/ex_data.md)
  [`ex_comp`](https://dragonfly-science.github.io/ComPoM/reference/ex_data.md)
  [`ex_catch`](https://dragonfly-science.github.io/ComPoM/reference/ex_data.md)
  : Simulated example data
