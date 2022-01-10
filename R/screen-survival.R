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
#' individual infected i days ago has been diagnosed.
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

  if(length(na.omit(tslt[sub_pop])) == 0 ||
     length(na.omit(ever_hiv_test[sub_pop])) == 0)
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
