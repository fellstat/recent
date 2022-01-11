## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(rita)
data("assay_data")
head(assay_data[1:7])

## -----------------------------------------------------------------------------
inc <- rita_incidence(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    elite_cntr=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    frr = lag_avidity_frr()[1],
    assay_surv = lag_avidity_survival(2 * 365)
  )
 knitr::kable(inc)

## -----------------------------------------------------------------------------
inc <- rita_incidence(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    elite_cntr=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    frr = lag_avidity_frr()[1],
    assay_surv = lag_avidity_survival(2 * 365),
    test_history_population = "negative"
  )
knitr::kable(inc)

## -----------------------------------------------------------------------------
rep_weights <-  dplyr::select(assay_data, dplyr::contains("btwt"))
ri <- rita_bootstrap(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    elite_cntr=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    rep_weights = rep_weights,
    rep_weight_type = "JK1"
  )
knitr::kable(ri)

