# Create a spatial phylogenetic object

This function creates a `phylospatial` object. This is the core data
type in the phylospatial library, and is a required input to most other
functions in the package. The two essential components of a spatial
phylogenetic object are a phylogenetic tree and an community data set.

## Usage

``` r
phylospatial(
  comm,
  tree = NULL,
  spatial = NULL,
  data_type = c("auto", "probability", "binary", "abundance", "other"),
  clade_fun = NULL,
  build = TRUE,
  check = TRUE,
  area_tol = 0.01
)
```

## Arguments

- comm:

  Community data representing the distribution of terminal taxa across
  sites. Can be a matrix with a column per terminal and a row per site,
  a
  [SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with one layer per terminal, or a `sf` data with a column per
  terminal. Taxa whose names do not match between column/layer names in
  `comm` and tip labels in `tree` will be dropped with a warning (unless
  `build = FALSE`).

- tree:

  Phylogeny of class
  [phylo](https://rdrr.io/pkg/ape/man/read.tree.html). Terminals whose
  names do not match `comm` will be dropped with a warning (unless
  `build = FALSE`). If this argument is not provided, terminals are
  assumed to follow a "star" tree with uniform branch lengths, which
  will lead to non-phylogenetic versions of any analyses done with the
  resulting `phylospatial` object. Must be a rooted tree.

- spatial:

  An optional `SpatRaster` layer or `sf` object indicating site
  locations. The number of cells or rows must match `comm`. Ignored if
  `comm` is a `SpatRaster` or `sf` object.

- data_type:

  Character giving the data type of `comm`. Must be "binary",
  "probability", "abundance", "auto" (the default), or "other". This
  determines how community values for clades are calculated from the
  values for terminal taxa. If "binary" (presence-absence), a clade is
  considered present in a site if any terminal in the clade is present.
  If "probability," clade probabilities are calculated as the
  probability that at least one terminal is present in a site. If
  "abundance," clade abundances are calculated as the sum of abundances
  for terminals in the clade in a site. If "auto," an attempt is made to
  guess which of these three data types was provided. This argument is
  ignored if `clade_fun` is provided, or if `build = FALSE`. If "other",
  a custom `clade_fun` must be supplied.

- clade_fun:

  Function to calculate the local community weight for a clade based on
  community weights for tips found in a given location. Must be either
  NULL (the default, in which case the default function for the selected
  `data_type` is used) or a summary function that takes a numeric vector
  and returns a single numeric output. Ignored if `comm` already
  includes clade ranges.

- build:

  Logical indicating whether `comm` already includes clade ranges that
  should be used instead of building new ones. Default is `TRUE`. If
  `FALSE`, `clade_fun` is ignored, no checks are performed to harmonize
  the tip labels and the community data, and the columns of `comm` must
  exactly match the order of `tree` edges including tips and larger
  clades. If clade ranges are included in `comm` but `build = TRUE`,
  they will be dropped and new clade ranges will be built.

- check:

  Logical indicating whether community data should be validated. Default
  is TRUE.

- area_tol:

  Numeric value giving tolerance for variation in the area of sites.
  Default is `0.01`. If the coefficient of variation in the area or
  length of spatial units (e.g. grid cells) exceeds this value, an error
  will result. This check is performed because various other functions
  in the library assume that sites are equal area. This argument is
  ignored if `check = FALSE` or if no spatial data is provided.

## Value

A `phylospatial` object, which is a list containing the following
elements:

- "data_type"::

  Character indicating the community data type

- "tree"::

  Phylogeny of class `phylo`

- "comm"::

  Community matrix, including a column for every terminal taxon and
  every larger clade. Column order corresponds to tree edge order.

- "spatial"::

  A `SpatRaster` or `sf` providing spatial coordinates for the rows in
  `comm`. May be missing if no spatial data was supplied.

- "dissim"::

  A community dissimilary matrix of class `dist` indicating pairwise
  phylogenetic dissimilarity between sites. Missing unless
  `ps_dissim(..., add = TRUE)` is called.

## Details

This function formats the input data as a `phylospatial` object. Beyond
validating, cleaning, and restructing the data, the main operation it
performs is to compute community occurrence data for every internal
clade on the tree. For a given clade and site, community data for all
the terminals in the clade are used to calculate the clade's occurrence
value in the site. As described above, this calculation can happen in
various ways, depending on what type of community data you have (e.g.
binary, probability, or abundance) and how you want to summarize them.
By default, the function tries to detect your `data_type` and use it to
automatically select an appropriate summary function as described above,
but you can override this by providing your own function to `clade_fun`.

You can also disable construction of the clade community matrix columns
altogether by setting `build = FALSE`). This is atypical, but you might
want to use this option if you have your own distribution data data on
all clades (e.g. from modeling occurrence probabilities for clades in
addition to terminal species), or if your community data comes from a
previously-constructed `phylospatial` object.

## Examples

``` r
# \donttest{
# load species distribution data and phylogeny
comm <- terra::rast(system.file("extdata", "moss_comm.tif", package = "phylospatial"))
tree <- ape::read.tree(system.file("extdata", "moss_tree.nex", package = "phylospatial"))

# construct `phylospatial` object
ps <- phylospatial(comm, tree)
#> Community data type detected: probability
ps
#> `phylospatial` object
#>   - 884 lineages across 1116 sites
#>   - community data type: probability 
#>   - spatial data class: SpatRaster 
#>   - dissimilarity data: none 
# }
```
