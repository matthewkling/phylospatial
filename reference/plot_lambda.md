# Plot alternative lambda values

Show a plot illustrating alternative values for the `lambda` parameter
in
[ps_prioritize](https://matthewkling.github.io/phylospatial/reference/ps_prioritize.md).
Lambda determines the shape of the "benefit" function that determines
the conservation value of protecting a given proportion of the
geographic range of a species or clade. Positive values place a higher
priority on protecting additional populations of largely unprotected
taxa, whereas negative values place a higher priority on protecting
additional populations of relatively well-protected taxa. The default
value used by
[ps_prioritize](https://matthewkling.github.io/phylospatial/reference/ps_prioritize.md)
is 1.

## Usage

``` r
plot_lambda(lambda = c(-1, -0.5, 0, 0.5, 2, 1))
```

## Arguments

- lambda:

  A vector of lambda values to plot

## Value

Plots a figure

## Examples

``` r
plot_lambda()

plot_lambda(seq(0, 3, .1))

```
