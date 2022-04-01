

#' Assay Based Incidence Estimation
#' @param recent Logical. Tests recent on assay.
#' @param undiagnosed Logical. No previous diagnosis.
#' @param low_viral Logical. Has low viral load (< 1000).
#' @param hiv Logical. Is HIV positive.
#' @param weights Survey weights.
#' @param tslt Time since last HIV test (days).
#' @param ever_hiv_test Subject has been tested for HIV in the past.
#' @param tau long term cut-off (years).
#' @param frr Reference false recency rate among treatment naive non-elite controller non-AIDS individuals.
#' @param test_history_population If undiagnosed, the testing histories of undiagnosed HIV+ people are used. If negative, the HIV- population is used.
#' @param assay_surv Survival function vector for assay among treatment naive non-elite controller non-AIDS individuals.
#' @param diag_surv time to diagnosis survival function vector. If specified, overrides the internal calculation.
#' @param treated A logical vector indicating a subject is on treatment. Only needed in the case of the use of RITA2 screening.
#' @param treat_surv Probability an individual diagnosed i days ago is not on treatment.
#' @details
#' This function estimates HIV incidence for cross-sectional survey designs using a recency assay
#' combined with a Recent Infection Testing Algorithm (RITA) screening step, which is used to
#' remove long-term individuals with elevated false recency rates on the assay. Two RITA algorithms
#' are supported.
#' RITA3 treats all individuals with either a previous diagnosis (as determined by self report or ARV biomarkers) or
#' a viral load <1,000 c/ml as non-recent regardless of the result of the recency assay. RITA2
#' treats all individuals who are either on treatment or have a viral load <1,000 c/ml as
#' non-recent.
#' The default RITA is RITA3. If 'treated' is non-null, RITA2 will be used. RITA2 also requires
#' a vector 'treat_surv' whose ith element represents the probability that an individual
#' diagnosed i days ago is not on treatment.
#' @returns
#' A data.frame with the following values:
#'
#' 1. `incidence`: The incidence.
#' 2. `residual_frr`: The false recency rate accounting for the screening process.
#' 3. `omega_rs`: The mean duration of recency up to tau accounting for the screening process.
#' 4. `P(R|S)` : The proportion of screened in individual who test recent.
#' 5. `P(S|H)` : The proportion of HIV+ individuals that are screened in.
#' 6. `P(H)` : HIV prevalence.
#'
#' @examples
#' data("assay_data")
#' rita_incidence(
#' recent=assay_data$recent,
#' undiagnosed=assay_data$undiagnosed,
#' low_viral=assay_data$elite_cntr,
#' hiv=assay_data$hiv,
#' weights=assay_data$weights,
#' tslt=assay_data$tslt,
#' ever_hiv_test=assay_data$ever_hiv_test
#' )
#'
#' # RITA2 Screening
#' ## Posit an average time from diagnosis to treatment of 150 days
#' treat_surv <- 1 - pexp(1:(365*2), 1/150)
#'
#' ## Create a dummy variable for treatment
#' assay_data$treated <- !assay_data$undiagnosed
#' assay_data$treated[assay_data$undiagnosed][c(40L, 47L, 59L, 63L, 83L, 157L, 164L, 166L, 194L, 209L)] <- FALSE
#'
#' # Calculate incidence using RITA2 screening (i.e. screen as non-recent if either treated or low viral load)
#' rita_incidence(
#' recent=assay_data$recent,
#' undiagnosed=assay_data$undiagnosed,
#' low_viral=assay_data$elite_cntr,
#' hiv=assay_data$hiv,
#' weights=assay_data$weights,
#' tslt=assay_data$tslt,
#' ever_hiv_test=assay_data$ever_hiv_test,
#' treated = assay_data$treated,
#' treat_surv = treat_surv
#' )
#' @export
rita_incidence <- function(
  recent,
  undiagnosed,
  low_viral,
  hiv,
  tslt,
  ever_hiv_test,
  weights = rep(1, length(recent)),
  tau = 2,
  frr = lag_avidity_frr()[1],
  test_history_population = c("undiagnosed", "negative"),
  assay_surv = lag_avidity_survival(tau * 365),
  diag_surv = NULL,
  treated = NULL,
  treat_surv = NULL
){
  test_history_population <- match.arg(test_history_population)
  tau_days <- tau * 365

  if(is.null(diag_surv))
    diag_surv <- diagnosis_survival(
      undiagnosed,
      tslt,
      ever_hiv_test,
      hiv,
      weights,
      n=tau_days,
      population = test_history_population
    )

  frr <- as.vector(frr)

  omega_s <- sum(diag_surv[1:tau_days]) / 365
  omega <- sum(assay_surv[1:tau_days] * diag_surv[1:tau_days]) / 365

  screen_in <- undiagnosed & (!low_viral)
  screen_in[!hiv] <- FALSE

  if(!is.null(treated)){
    if(length(treated) != length(undiagnosed))
      stop("length(treated) != length(undiagnosed)")
    if(any(treated[undiagnosed]))
      stop("Undiagnosed individuals cannot be on treatment")
    diag_treat_surv <- apply_time_to_treatment(diag_surv, treat_surv)
    omega_s <- sum(diag_treat_surv[1:tau_days]) / 365
    omega <- sum(assay_surv[1:tau_days] * diag_treat_surv[1:tau_days]) / 365
    screen_in <- (!treated) & (!low_viral)
    screen_in[!hiv] <- FALSE
  }

  phiv <- sum(hiv * weights, na.rm=TRUE) /
      sum(weights[!is.na(hiv)])


  pscreen <- sum(screen_in * weights, na.rm=TRUE) /
      sum(weights[!is.na(screen_in)])

  #precent_ghiv <- sum(recent * screen_in * hiv * weights, na.rm=TRUE) /
  #    sum((!is.na(recent * screen_in * hiv)) * hiv * weights)

  precent_gscreen <- sum(recent * screen_in * weights, na.rm=TRUE) /
      sum( (screen_in * weights)[!is.na(recent * screen_in)])

  lambda <- (precent_gscreen - frr) * pscreen / ((1-phiv) * (omega - frr * omega_s))
  rita_frr <- frr * (pscreen - omega_s * lambda * (1-phiv)) /
    (phiv - tau  * lambda * (1-phiv))

  data.frame(
    incidence = lambda,
    residual_frr = rita_frr,
    omega_rs = omega,
    omega_s = omega_s,
    `P(R|S)` = precent_gscreen,
    `P(S|H)` = pscreen / phiv,
    `P(H)` = phiv,
    check.names=FALSE
  )
}


#' Survey bootstrap
#' @param recent Logical. Tests recent on assay.
#' @param undiagnosed Logical. No previous diagnosis.
#' @param low_viral Logical.  Has low viral load (< 1000).
#' @param hiv Logical. Is HIV positive.
#' @param weights Survey weights.
#' @param tslt Time since last HIV test (days).
#' @param ever_hiv_test Subject has been tested for HIV in the past.
#' @param tau long term cut-off (years).
#' @param frr False recency rate among treatment naive non-elite controller non-AIDS individuals.
#' @param test_history_population If undiagnosed, the testing histories of undiagnosed HIV+ people are used. If negative, the HIV- population is used.
#' @param assay_surv Survival function vector for assay among treatment naive non-elite controller non-AIDS individuals.
#' @param diag_surv time to diagnosis survival function vector.
#' @param treated A logical vector indicating a subject is on treatment. Only needed in the case of the use of RITA2 screening.
#' @param treat_surv Probability an individual diagnosed i days ago is not on treatment.
#' @param rep_weights A data.frame of replicate weights. See survey::svrrepdesign
#' @param rep_weight_type The type of resampling weights. See svrepdesign.
#' @param combined_weights TRUE if the rep_weights already include the sampling weights. This is usually the case.
#' @param conf_level confidence level for bootstrap interval.
#' @param show_progress If TRUE, prints bootstrap progress. This may also be a callback function taking one parameter equal to the index of the current replicate.
#' @param ... additional parameters to svrepdesign.
#' @returns
#' A data.frame with columns for the estimate, standard error, lower confidence bound and upper confidence bound.
#' Rows are defined by:
#'
#' 1. `incidence`: The incidence.
#' 2. `residual_frr`: The false recency rate accounting for the screening process.
#' 3. `omega_rs`: The mean duration of recency up to tau accounting for the screening process.
#' 4. `P(R|S)` : The proportion of screened in individual who test recent.
#' 5. `P(S|H)` : The proportion of HIV+ individuals that are screened in.
#' 6. `P(H)` : HIV prevalence.
#' @examples
#' data("assay_data")
#' rep_weights <-  dplyr::select(assay_data, dplyr::contains("btwt"))
#' rita_bootstrap(
#' recent=assay_data$recent,
#' undiagnosed=assay_data$undiagnosed,
#' low_viral=assay_data$elite_cntr,
#' hiv=assay_data$hiv,
#' weights=assay_data$weights,
#' tslt=assay_data$tslt,
#' ever_hiv_test=assay_data$ever_hiv_test,
#' rep_weights = rep_weights,
#' rep_weight_type = "JK1"
#' )
#' @export
rita_bootstrap <- function(
  recent,
  undiagnosed,
  low_viral,
  hiv,
  tslt,
  ever_hiv_test,
  weights,
  rep_weights = NULL,
  rep_weight_type = c("BRR", "Fay", "JK1", "JK2", "JKn", "bootstrap", "other"),
  combined_weights = TRUE,
  tau = 2,
  frr = lag_avidity_frr()[1],
  test_history_population = c("undiagnosed", "negative"),
  assay_surv = lag_avidity_survival(tau * 365),
  diag_surv = NULL,
  treated = NULL,
  treat_surv = NULL,
  conf_level=.95,
  show_progress = TRUE,
  ...
){
  test_history_population <- match.arg(test_history_population)
  fun <- function(wts){
    as.matrix(rita_incidence(
      recent=recent,
      undiagnosed = undiagnosed,
      low_viral = low_viral,
      hiv = hiv,
      tslt = tslt,
      ever_hiv_test = ever_hiv_test,
      weights = wts,
      tau = tau,
      frr = frr,
      test_history_population = test_history_population,
      assay_surv = assay_surv,
      diag_surv = diag_surv,
      treated = treated,
      treat_surv = treat_surv
    ))
  }
  values <- fun(weights)

  nr <- nrow(values)
  nc <- ncol(values)

  not_miss <- (rowSums(is.na(rep_weights)) + is.na(weights)) < 0.5
  if(sum(not_miss) < 2)
    stop("Too few observations with non-missing rep_weights and weights")
  recent <- recent[not_miss]
  undiagnosed <- undiagnosed[not_miss]
  low_viral <- low_viral[not_miss]
  hiv <- hiv[not_miss]
  tslt <-tslt[not_miss]
  ever_hiv_test <- ever_hiv_test[not_miss]
  weights <- weights[not_miss]
  rep_weights <- rep_weights[not_miss,]
  if(!is.null(treated))
    treated <- treated[not_miss]

  rep_design <- survey::svrepdesign(repweights = rep_weights,
                                    weights = weights,
                                    data = rep_weights,
                                    combined.weights = combined_weights,
                                    type = rep_weight_type,
                                    ...)
  scale <- rep_design$scale
  rscales <- rep_design$rscales
  mse <- rep_design$mse
  rep_weights <- weights(rep_design)
  nrep <- ncol(rep_weights)

  errors <- list()
  estimates <- array(NA, dim=c(nr, nc, nrep))
  if(!is.function(show_progress) && show_progress)
    prog_bar <- progress::progress_bar$new(total=nrep)
  for(i in 1:nrep){
    if(is.function(show_progress))
      show_progress(i)
    else if(show_progress)
      prog_bar$tick()
    val <- try(fun(rep_weights[,i]))
    if(!inherits(val, "try-error"))
      estimates[,,i] <- val#fun(rep_weights[,i])
    else
      errors[[length(errors) + 1]] <- val
  }

  vars <- matrix(NA,nrow=nr,ncol=nc)
  for(i in 1:nr){
    for(j in 1:nc){
      if(!all(is.na(estimates[i,j,])))
        vars[i,j] <- survey::svrVar(
          estimates[i,j,],
          scale,
          rscales,
          mse = mse,
          coef = values[i,j]
          )
    }
  }
  bb <- list(value=values,
             var = vars,
             replicates = estimates,
             nrep = nrep,
             errors = errors)
  res <- .table_rita_boot_est(bb, conf_level = conf_level)
  attr(res,"bootstraps") <- bb
  res
}





.table_rita_boot_est <- function(object, conf_level=.95, ...){
  cival <- stats::qnorm(1 - (1 - conf_level) / 2)
  res <- as.data.frame(t(rbind(
    object$value,
    sqrt(object$var),
    object$value - cival * sqrt(object$var),
    object$value + cival * sqrt(object$var)
  )))
  colnames(res) <- c(
    "estimate",
    "std_error",
    "lower_bound",
    "upper_bound"
  )
  res
}

