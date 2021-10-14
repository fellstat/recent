

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
#' @param assay_surv Survival function vector for assay among treatment naive non-elite controller non-AIDS individuals
#' @param aids_surv Survival function vector for time to AIDS in untreated individuals.
#' @export
rita_incidence <- function(
  recent,
  undiagnosed,
  elite_cntr,
  hiv,
  tslt,
  ever_hiv_test,
  weights = rep(1, length(recent)),
  diag_surv = NULL,
  tau = 2,
  frr = lag_avidity_frr()[1],
  assay_surv = lag_avidity_survival(tau * 365),
  aids_surv = aids_survival(tau * 365)
){
  tau_days <- tau * 365
  if(is.null(diag_surv)){
    if(is.null(tslt))
      stop("tslt is required if diag_surv is not specified")
    if(is.null(ever_hiv_test))
      stop("ever_hiv_test is required if diag_surv is not specified")
    diag_surv <- diagnosis_survival(undiagnosed, tslt, ever_hiv_test, weights, n=tau_days)
  }
  frr <- as.vector(frr)
  aid_surv <- aids_survival(tau_days)
  omega_s <- sum(diag_surv[1:tau_days] * aids_surv[1:tau_days]) / 365

  omega <- sum(assay_surv[1:tau_days] * diag_surv[1:tau_days] * aids_surv[1:tau_days]) / 365

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
    (phiv - tau  * lambda * (1-phiv))

  data.frame(
    lambda = lambda,
    rita_frr = rita_frr,
    omega = omega ,
    omega_s = omega_s,
    `P(R|S)` = precent_gscreen,
    `P(S|H)` = pscreen / phiv,
    `P(H)` = phiv,
    check.names=FALSE
  )
}


#' Survey bootstrap
#' @export
rita_bootstrap <- function(
  recent,
  undiagnosed,
  elite_cntr,
  hiv,
  tslt,
  ever_hiv_test,
  weights,
  rep_weights = NULL,
  rep_weight_type=c("BRR", "Fay", "JK1","JK2","JKn","bootstrap","other"),
  combined_weights=TRUE,
  diag_surv = NULL,
  tau = 2,
  frr = lag_avidity_frr()[1],
  assay_surv = lag_avidity_survival(tau * 365),
  aids_surv = aids_survival(tau * 365),
  show_progress = TRUE,
  ...
){
  fun <- function(wts){
    as.matrix(rita_incidence(
      recent=recent,
      undiagnosed = undiagnosed,
      elite_cntr = elite_cntr,
      hiv = hiv,
      tslt=tslt,
      ever_hiv_test = ever_hiv_test,
      weights = wts,
      diag_surv = diag_surv,
      tau = tau,
      frr = frr,
      assay_surv = assay_surv,
      aids_surv = aids_surv
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
  elite_cntr <- elite_cntr[not_miss]
  hiv <- hiv[not_miss]
  tslt <-tslt[not_miss]
  ever_hiv_test <- ever_hiv_test[not_miss]
  weights <- weights[not_miss]
  rep_weights <- rep_weights[not_miss,]

  rep_design <- survey::svrepdesign(repweights=rep_weights,weights=weights,
                                    data=rep_weights,
                                    combined.weights=combined_weights,
                                    type=rep_weight_type,
                                    ...)
  scale <- rep_design$scale
  rscales <- rep_design$rscales
  mse <- rep_design$mse

  rep_weights <- weights(rep_design)
  nrep <- ncol(rep_weights)

  errors <- list()
  estimates <- array(NA, dim=c(nr, nc, nrep))
  for(i in 1:nrep){
    if(is.function(show_progress))
      show_progress(i)
    else if(show_progress)
      cat(".")
    val <- try(fun(rep_weights[,i]))
    if(!inherits(val, "try-error"))
      estimates[,,i] <- val#fun(rep_weights[,i])
    else
      errors[[length(errors) + 1]] <- val
  }
  if(!is.function(show_progress) && show_progress)
    cat("\n")
  vars <- matrix(NA,nrow=nr,ncol=nc)
  for(i in 1:nr){
    for(j in 1:nc){
      if(!all(is.na(estimates[i,j,])))
        vars[i,j] <- survey::svrVar(estimates[i,j,], scale, rscales, mse=mse,coef=values[i,j])
    }
  }
  bb <- list(value=values,
             var = vars,
             replicates=estimates,
             nrep=nrep,
             errors=errors)
  class(bb) <- c("rita_boot_est","list")
  bb
}





#' print bootstrap
#' @param x a inc_boot_est object
#' @param ... additional parameters for summary.inc_boot_est
#' @export
print.rita_boot_est <- function(x, ...){
  res <- summary(x,...)
  print(res)
}

#' summary of incidence bootstrap
#' @param object a inc_boot_est object
#' @param conf_level confidence level for intervals
#' @param ... additional parameters for print.data.frame
#' @export
summary.rita_boot_est <- function(object, conf_level=.95, ...){
  cival <- qnorm(1 - (1 - conf_level) / 2)
  res <- as.data.frame(t(rbind(
    object$value,
    sqrt(object$var),
    object$value - cival * sqrt(object$var),
    object$value + cival * sqrt(object$var)
  )))
  colnames(res) <- c(
    "estimate",
    "se",
    paste0(100*conf_level," CI lower"),
    paste0(100*conf_level," CI upper")
  )
  res
}
