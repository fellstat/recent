
# rita: Incidence Estimation Using Recency Assays

<!-- badges: start -->
<!-- badges: end -->

Provides functions and methods to estimate incidence from a cross-sectional survey. A Recent Infection Testing Algorithm (RITA) is used to identify likely recent infections and the methods of Fellows (2021) are used to calculate incidence based on the frequency of these likely recent infections.

## Installation

You can install from github using

``` r
devtools::install_github("fellstat/rita",auth_token="<your github auth token>")
```

# Usage

**[A brief tour of the package](https://htmlpreview.github.io/?https://github.com/fellstat/rita/blob/main/inst/doc/introduction.html)**

## Example

``` r
library(rita)
data(assay_recent)
rita_incidence(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    elite_cntr=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test
  )
```

