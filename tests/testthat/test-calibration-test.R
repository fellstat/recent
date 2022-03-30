test_that("diagnosis_survival", {
  data("assay_data")
  ds <- diagnosis_survival(
    assay_data$undiagnosed,
    assay_data$tslt,
    assay_data$ever_hiv_test,
    assay_data$hiv,
    assay_data$weights,
    n=40)
  ds_ref <- c(0.999999999254463, 0.999999995735573, 0.999999988172136, 0.999999975607732,
              0.999999957236054, 0.999999932345329, 0.985622917854146, 0.964535215673545,
              0.964535169254846, 0.964535114315464, 0.964535050378901, 0.964534976991274,
              0.96453489371831, 0.964534800142956, 0.964534695863445, 0.964534580491691,
              0.964534453651958, 0.964534314979725, 0.964534164120719, 0.964534000730082,
              0.964533824471638, 0.964533635017262, 0.964533432046313, 0.964533215245133,
              0.964532984306605, 0.957876340289432, 0.946198913238991, 0.936192994863821,
              0.932983387017365, 0.932983090997963, 0.932982779634702, 0.912214407934872,
              0.906356803565437, 0.896104965132648, 0.878844048302637, 0.87435042272391,
              0.870817150394265, 0.870816752884103, 0.861533929404766, 0.861533503901948
  )
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
  ds2_ref <- c(0.999999999254463, 0.999999995735573, 0.999999988172136, 0.999999975607732,
               0.999999957236054, 0.997393491868288, 0.980695341584368, 0.961521881541127,
               0.945954844457079, 0.945954790576022, 0.945954727871102, 0.945954655897177,
               0.945954574228342, 0.945954482455577, 0.945954380184854, 0.945954267035565,
               0.945954142639209, 0.945954006638289, 0.945953858685356, 0.945953698442193,
               0.945953525579105, 0.945953339774284, 0.945953140713266, 0.944867333052843,
               0.943782365530925, 0.940736990723371, 0.934281398312705, 0.926027197790439,
               0.913849902873736, 0.901793002582042, 0.892154891263378, 0.884751453644247,
               0.871908202441554, 0.863604313099702, 0.855382578669327, 0.846952861813802,
               0.835338331532795, 0.827749295970403, 0.823676276276591, 0.823675869471281
  )
  expect_true(sum(abs(ds2 - ds2_ref))<.000001)
})



test_that("rita_incidence", {
  data("assay_data")
  ri <- rita_incidence(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    low_viral=assay_data$elite_cntr,
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
    low_viral=assay_data$elite_cntr,
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
    low_viral=assay_data$elite_cntr,
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
    low_viral=assay_data$elite_cntr,
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



test_that("RITA2", {
  data("assay_data")
  diag_surv <- diagnosis_survival(
    assay_data$undiagnosed,
    assay_data$tslt,
    assay_data$ever_hiv_test,
    assay_data$hiv,
    assay_data$weights,
    n=365*2)

  # Posit an average time to treatment of 150 days
  treat_surv <- 1 - pexp(1:(365*2), 1/150)

  # Calculate the treatment survival function
  diag_treat_surv <- apply_time_to_treatment(diag_surv, treat_surv)

  # Compare survival curve for time to diagnosis (red) vs time to treatment (black)
  plot(diag_treat_surv, type="l")
  points(diag_surv, type="l",col="red")

  #Create a dummy variable for treatment
  assay_data$treated <- assay_data$undiagnosed
  assay_data$treated[assay_data$undiagnosed][c(40L, 47L, 59L, 63L, 83L, 157L, 164L, 166L, 194L, 209L)] <- FALSE

  #Calculate incidence
  ri1 <- rita_incidence(
    recent=assay_data$recent,
    undiagnosed=!assay_data$treated, #used treated indicator in place of undiagnosed for screening
    low_viral=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    diag_surv = diag_treat_surv
  )

  ri2 <- rita_incidence(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    low_viral=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    treated = assay_data$treated,
    treat_surv = treat_surv
  )
  expect_identical(ri1,ri2)

  ri1_ref <- structure(list(incidence = 0.00383206717870841, residual_frr = 0.000842313332886876,
                            omega_rs = 0.335938952038046, omega_s = 1.33546256766089,
                            `P(R|S)` = 0.0287379226230863, `P(S|H)` = 0.164430714369552,
                            `P(H)` = 0.24804274938758), class = "data.frame", row.names = c(NA,
                                                                                            -1L))
  expect_true(sum(abs(ri1 - ri1_ref))<.000001)

  expect_warning(rboot <- rita_bootstrap(
    recent=assay_data$recent,
    undiagnosed=assay_data$undiagnosed,
    low_viral=assay_data$elite_cntr,
    hiv=assay_data$hiv,
    weights=assay_data$weights,
    tslt=assay_data$tslt,
    ever_hiv_test=assay_data$ever_hiv_test,
    rep_weights = assay_data %>% dplyr::select(contains("btwt")),
    rep_weight_type = "JK2",
    treated = assay_data$treated,
    treat_surv = treat_surv,
    show_progress = FALSE
  ))
  expect_true(sum(abs(ri1 - rboot[,1]))<.000001)

})
