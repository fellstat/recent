#' Survival probability for RITA Screening
#' @param undiagnosed A logical vector indicating if an individual has never been diagnosed as HIV+.
#' @param tslt Time since last test in days
#' @param ever_hiv_test A logical vector indicating whether an individual has ever had an HIV test
#' @param hiv HIV positive status.
#' @param weights Survey weights
#' @param n Number of days to calculate
#' @param population If undiagnosed, the testing histories of undiagnosed HIV+ people are used. If negative, the HIV- population is used.
#' @returns
#' A vector where the ith element is the probability that an
#' individual infected i days ago has not been diagnosed.
#' @examples
#'   data("assay_data")
#'   diagnosis_survival(
#'     assay_data$undiagnosed,
#'     assay_data$tslt,
#'     assay_data$ever_hiv_test,
#'     assay_data$hiv,
#'     assay_data$weights,
#'     n=40)
#' @export
diagnosis_survival <- function(
  undiagnosed,
  tslt,
  ever_hiv_test,
  hiv,
  weights,
  n = 365*5,
  population = c("undiagnosed", "negative")){

  population <- match.arg(population)
  if(population == "undiagnosed"){
    sub_pop <- undiagnosed
    sub_pop[!hiv] <- FALSE
  }else{
    sub_pop <- !hiv
  }
  sub_pop[is.na(sub_pop)] <- FALSE

  if(length(stats::na.omit(tslt[sub_pop])) == 0 ||
     length(stats::na.omit(ever_hiv_test[sub_pop])) == 0)
    stop("diagnosis_survival: no non-missing testing histories")

  tst_surv <- sapply(1:n, function(x){
    sum((tslt[sub_pop] > x) * weights[sub_pop], na.rm=TRUE) /
      sum((!is.na(tslt[sub_pop])) * weights[sub_pop])
  })
  p_never_test <- sum((!ever_hiv_test[sub_pop]) * weights[sub_pop], na.rm=TRUE) /
    sum((!is.na(ever_hiv_test[sub_pop])) * weights[sub_pop])
  tst_surv <- 1 - (1 - tst_surv) * (1 - p_never_test)
  tst_surv
}

#' Obtain Probability Curve for Treatment
#' @param diag_surv Time to diagnosis probability function vector. The probability an individual diagnosed i days ago is not diagnosed.
#' @param treat_surv Time to treatment. The probability an individual diagnosed i days ago is not on treatment.
#' @details
#' If the RITA algorithm excludes individuals based on treatment rather than diagnosis, this function can be used to calculate the survival function for use in icidence estimation.
#' @returns
#' A vector where the ith element is the probability that an
#' individual infected i days ago has either not been diagnosed or is not on treatment.
#' @examples
#' data("assay_data")
#'
#' #Obtain the diagnosis survival function
#' diag_surv <- diagnosis_survival(
#'   assay_data$undiagnosed,
#'   assay_data$tslt,
#'   assay_data$ever_hiv_test,
#'   assay_data$hiv,
#'   assay_data$weights,
#'   n=365*2)
#'
#' # Posit an average time to treatment of 150 days
#' treat_surv <- 1 - pexp(1:(365*2), 1/150)
#'
#' # Calculate the treatment survival function
#' diag_treat_surv <- apply_time_to_treatment(diag_surv, treat_surv)
#'
#' # Compare survival curve for time to diagnosis (red) vs time to treatment (black)
#' plot(diag_treat_surv, type="l")
#' points(diag_surv, type="l",col="red")
#'
#' #Create a dummy variable for treatment
#' assay_data$treated <- assay_data$undiagnosed
#'
#' #Calculate incidence
#' rita_incidence(
#' recent=assay_data$recent,
#' undiagnosed=assay_data$treated, #used treated indicator in place of undiagnosed for screening
#' elite_cntr=assay_data$elite_cntr,
#' hiv=assay_data$hiv,
#' weights=assay_data$weights,
#' tslt=assay_data$tslt,
#' ever_hiv_test=assay_data$ever_hiv_test,
#' diag_surv = diag_treat_surv
#' )
#' @export
apply_time_to_treatment <- function(diag_surv, treat_surv){
  n <- length(diag_surv)
  if(n < 1)
    stop("length(diag_surv) < 1")
  if(length(treat_surv) < n)
    stop("The length of treat_surv must be at least as long as diag_surv")
  diag_haz <- diff(diag_surv)
  res <- rep(0, n)
  for(i in 1:(n-1)){
    res[i:n] <- res[i:n] + diag_haz[i] * (1-treat_surv[1:(n - i + 1)])
  }
  1 + res
}
