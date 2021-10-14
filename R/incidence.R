

#' Assay Based Incidence Estimation
#' @param recent Logical. Tests recent on assay
#' @param undiagnosed Logical. No previous diagnosis
#' @param elite_cntr Logical. Is an elite controller (VL < 1000)
#' @param hiv Logical. Is HIV positive
#' @param weights Survey weights
#' @param tslt Time since last HIV test
#' @param ever_hiv_test Subject has been tested for HIV in the past
#' @param diag_surv time to diagnosis survival function vector
#' @param tau long term cut-off
#' @param frr False recency rate among treatment naive non-elite controller non-AIDS individuals
#' @param frr_se Standard error of FRR
#' @param assay_surv Survival function vector for assay among treatment naive non-elite controller non-AIDS individuals
#' @param aids_surv Survival function vector for time to AIDS in untreated individuals.
#' @export
rita_incidence <- function(
  recent,
  undiagnosed,
  elite_cntr,
  hiv,
  weights,
  tslt = NULL,
  ever_hiv_test = NULL,
  diag_surv = NULL,
  tau = 2 * 365,
  frr = lag_avidity_frr()[1],
  frr_se = lag_avidity_frr()[2],
  assay_surv = lag_avidity_survival(),
  aids_surv = aids_survival(tau)
){
  if(is.null(diag_surv)){
    if(is.null(tslt))
      stop("tslt is required if diag_surv is not specified")
    if(is.null(ever_hiv_test))
      stop("ever_hiv_test is required if diag_surv is not specified")
    diag_surv <- diagnosis_survival(undiagnosed, tslt, ever_hiv_test, weights, n=tau)
  }
  frr <- as.vector(frr)
  aid_surv <- aids_survival(tau)
  omega_s <- sum(diag_surv[1:tau] * aids_surv[1:tau]) / 365

  omega <- sum(assay_surv[1:tau] * diag_surv[1:tau] * aids_surv[1:tau]) / 365

  screen_in <- undiagnosed * (!elite_cntr)

  phiv <- sum(hiv * weights, na.rm=TRUE) /
      sum((!is.na(hiv)) * weights)


  pscreen <- sum(screen_in * weights, na.rm=TRUE) /
      sum((!is.na(screen_in)) * weights)

  #precent_ghiv <- sum(recent * screen_in * hiv * weights, na.rm=TRUE) /
  #    sum((!is.na(recent * screen_in * hiv)) * hiv * weights)

  precent_gscreen <- sum(recent * screen_in * weights, na.rm=TRUE) /
      sum((!is.na(recent * screen_in)) * screen_in * weights)

  lambda <- (precent_gscreen - frr) * pscreen / ((1-phiv) * (omega - frr * omega_s))
  rita_frr <- frr * (pscreen - omega_s * lambda * (1-phiv)) /
    (phiv - (tau / 365) * lambda * (1-phiv))

  data.frame(
    lambda = lambda,
    rita_frr = rita_frr,
    omega = omega * 365,
    omega_s = omega_s * 365
  )
}
