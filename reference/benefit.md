# Calculate taxon conservation benefit

Nonlinear function that converts proportion of range conserved into
conservation "benefit."

## Usage

``` r
benefit(x, lambda = 1)
```

## Arguments

- x:

  Fraction of taxon range protected (value between 0 and 1).

- lambda:

  Shape parameter.

## Value

Value between 0 and 1.
