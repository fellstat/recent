#' Survival probability for RITA Screening
#' @param undiagnosed A logical vector indicating if an individual has never been diagnosed as HIV+.
#' @param tslt Time since last test in days
#' @param ever_hiv_test A logical vector indicating whether an individual has ever had an HIV test
#' @param weights Survey weights
#' @param n Number of days to calculate
#' @returns
#' A vector where the ith element is the probability that an
#' individual infected i days ago has been diagnosed.
#' @export
diagnosis_survival <- function(undiagnosed, tslt, ever_hiv_test, weights, n=365*5){
  tst_surv <- sapply(1:n, function(x){
    sum((tslt[undiagnosed] > x) * weights[undiagnosed], na.rm=TRUE) /
      sum((!is.na(tslt[undiagnosed])) * weights[undiagnosed])
  })
  p_never_test <- sum((!ever_hiv_test[undiagnosed]) * weights[undiagnosed], na.rm=TRUE) /
    sum((!is.na(ever_hiv_test[undiagnosed])) * weights[undiagnosed])
  tst_surv <- 1 - (1 - tst_surv) * (1 - p_never_test)
  tst_surv
}
