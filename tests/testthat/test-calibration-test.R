test_that("diagnosis_survival", {
  data("assay_data")
  ds <- diagnosis_survival(
    assay_data$undiagnosed,
    assay_data$tslt,
    assay_data$ever_hiv_test,
    assay_data$hiv,
    assay_data$weights,
    n=40)
  ds_ref <- c(1, 1, 1, 1, 1, 1, 0.98562301612998, 0.964535350247737, 0.964535350247737,
              0.964535350247737, 0.964535350247737, 0.964535350247737, 0.964535350247737,
              0.964535350247737, 0.964535350247737, 0.964535350247737, 0.964535350247737,
              0.964535350247737, 0.964535350247737, 0.964535350247737, 0.964535350247737,
              0.964535350247737, 0.964535350247737, 0.964535350247737, 0.964535350247737,
              0.95787893358629, 0.946201730088066, 0.936196048975817, 0.932986711601312,
              0.932986711601312, 0.932986711601312, 0.912218572067184, 0.906361274008822,
              0.896109729774375, 0.878849074708043, 0.874355790732545, 0.870822878264839,
              0.870822878264839, 0.861540398772289, 0.861540398772289)
  expect_true(sum(abs(ds - ds_ref))<.000001)
  ds2 <- diagnosis_survival(
    assay_data$undiagnosed,
    assay_data$tslt,
    assay_data$ever_hiv_test,
    assay_data$hiv,
    assay_data$weights,
    n=40,
    population = "negative"
    )
  ds2_ref <- c(1, 1, 1, 1, 1, 0.997393559346621, 0.980695439368877, 0.961522015694892,
               0.945955021963414, 0.945955021963414, 0.945955021963414, 0.945955021963414,
               0.945955021963414, 0.945955021963414, 0.945955021963414, 0.945955021963414,
               0.945955021963414, 0.945955021963414, 0.945955021963414, 0.945955021963414,
               0.945955021963414, 0.945955021963414, 0.945955021963414, 0.944869424524847,
               0.943784680572044, 0.940739537618185, 0.93428417968315, 0.926030218738885,
               0.913853159277611, 0.901796502146821, 0.89215865116476, 0.88475549241198,
               0.871912502973354, 0.863608904933591, 0.855387470890678, 0.846958061617173,
               0.835343826038554, 0.827755118411892, 0.823682461366196, 0.823682461366196
  )
  expect_true(sum(abs(ds2 - ds2_ref))<.000001)
})



test_that("rita_incidence", {
  data("assay_data")
  ri <- rita_incidence(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    elite_cntr=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test
  )
  ri_ref <- structure(list(lambda = 0.0151573229739252, residual_frr = 0.000882682603705333,
                           omega_rs = 0.292134201319063, omega_s = 1.09313642218944, `P(R|S)` = 0.0728330452132411,
                           `P(S|H)` = 0.195313010675028, `P(H)` = 0.24804274938758), class = "data.frame", row.names = c(NA,
                                                                                                                         -1L))
  expect_true(sum(abs(ri - ri_ref))<.000001)

  dsurv <- diagnosis_survival(
    assay_data$undiagnosed,
    assay_data$tslt,
    assay_data$ever_hiv_test,
    assay_data$hiv,
    assay_data$weights,
    n=2*365)
  ri <- rita_incidence(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    elite_cntr=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    test_history_population = "negative",
    diag_surv = dsurv
  )
  expect_true(sum(abs(ri - ri_ref))<.000001)

  ri2 <- rita_incidence(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    elite_cntr=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    test_history_population = "negative"
  )
  ri2_ref <- structure(list(incidence = 0.0171628361526673, residual_frr = 0.000937839337768769,
                            omega_rs = 0.257254220714404, omega_s = 0.830829416650757,
                            `P(R|S)` = 0.0728330452132411, `P(S|H)` = 0.195313010675028,
                            `P(H)` = 0.24804274938758), class = "data.frame", row.names = c(NA,
             -1L))
  expect_true(sum(abs(ri2 - ri2_ref))<.000001)

})

test_that("rita_bootstrap", {
  data("assay_data")
  expect_warning(ri <- rita_bootstrap(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    elite_cntr=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    rep_weights = assay_data %>% dplyr::select(contains("btwt")),
    rep_weight_type = "JK2",
    show_progress = FALSE
  ))
#expect_true(sum(abs(ri - ri_ref))<.000001)
})
