test_that("diagnosis_survival", {
  data("assay_data")
  ds <- diagnosis_survival(assay_data$undiagnosed, assay_data$tslt, assay_data$ever_hiv_test, assay_data$weights,n=40)
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
  ri_ref <- structure(list(lambda = 0.0151573229739252, rita_frr = 0.000882682603705333,
                           omega = 0.292134201319063, omega_s = 1.09313642218944, `P(R|S)` = 0.0728330452132411,
                           `P(S|H)` = 0.195313010675028, `P(H)` = 0.24804274938758), class = "data.frame", row.names = c(NA,
                                                                                                                         -1L))
  expect_true(sum(abs(ri - ri_ref))<.000001)
})

test_that("rita_bootstrap", {
  data("assay_data")
  ri <- rita_bootstrap(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    elite_cntr=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    rep_weights = assay_data %>% dplyr::select(contains("btwt")),
    rep_weight_type = "JK2"
  )
                                                                                                     -1L))
  #expect_true(sum(abs(ri - ri_ref))<.000001)
})
