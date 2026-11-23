# ComPoM: Poisson-multinomial models to standardise and scale compositions

ComPoM fits Poisson-multinomial models to compositional count data (e.g.
length or age frequencies) and uses them to predict compositions for
strata with known catch. The predictions are aggregated into
catch-weighted compositions with uncertainty.

## Details

The workflow consists of four steps:

1.  [`data_prep()`](https://dragonfly-science.github.io/ComPoM/reference/data_prep.md):
    aggregate counts to strata x bins, fill zero counts and add stratum
    totals.

2.  [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md):
    fit the model with `brms` or `sdmTMB`, optionally with smooths of
    continuous covariates and bin-specific spatial and spatio-temporal
    random fields.

3.  Diagnostics and effects:
    [`post_pred_group()`](https://dragonfly-science.github.io/ComPoM/reference/post_pred_group.md),
    [`Fx_plot()`](https://dragonfly-science.github.io/ComPoM/reference/Fx_plot.md),
    [`cvar_smooth_plot()`](https://dragonfly-science.github.io/ComPoM/reference/cvar_smooth_plot.md),
    [`spatial_field_plot()`](https://dragonfly-science.github.io/ComPoM/reference/spatial_field_plot.md).

4.  Scaling:
    [`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md),
    then
    [`scaled_comp_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_comp_plot.md)
    or
    [`scaled_ridge_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_ridge_plot.md).

See
[`vignette("ComPoM")`](https://dragonfly-science.github.io/ComPoM/articles/ComPoM.md)
for the model equations and an example.

## See also

Useful links:

- <https://github.com/dragonfly-science/ComPoM>

- <https://dragonfly-science.github.io/ComPoM/>

- Report bugs at <https://github.com/dragonfly-science/ComPoM/issues>

## Author

**Maintainer**: Philipp Neubauer <philipp@dragonfly.co.nz> \[copyright
holder\]

Authors:

- Philipp Neubauer <philipp@dragonfly.co.nz> \[copyright holder\]
