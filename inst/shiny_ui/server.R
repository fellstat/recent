req <- function(pkg){
  r <- require(pkg, character.only = TRUE)
  if(!r){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
    return(2)
  }
  return(1)
}

req("shiny")
#req("promises")
#req("future")
# future::plan(future::multisession)
#req("ipc")
req("ggplot2")
req("rita")
req("survey")
req("shinyhelper")
options(shiny.maxRequestSize=300*1024^2)

shinyServer(function(input, output, session) {
  observe_helpers()

  get_raw_data <- reactive({
    inFile <- input$file1

    if (is.null(inFile))
      return(NULL)

    dat <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
    vars <- names(dat)

    for(name in c("hiv","recent","low_viral",
                  "undiagnosed", "ever_test","last_test","weights","strata","treated"))
      updateSelectizeInput(session, name, choices = c("Choose Variable"="",vars))
    updatePickerInput(session, "rep_weights", choices = vars)
    dat
  })

  get_numeric <- function(name){
    dat <- get_raw_data()
    if(is.null(dat) || input[[name]] == "")
      return(NULL)
    as.numeric(dat[[input[[name]]]])
  }

  get_categorical <- function(name){
    dat <- get_raw_data()
    if(is.null(dat) || input[[name]] == "")
      return(NULL)
    as.factor(dat[[input[[name]]]])
  }

  get_logical <- function(name){
    dat <- get_raw_data()
    if(is.null(dat) || input[[name]] == "")
      return(NULL)
    variable <- dat[[input[[name]]]]
    if(is.logical(variable))
      return(variable)
    if(is.numeric(variable))
      return(variable == max(variable, na.rm=TRUE))
    variable <- as.factor(variable)
    variable == max(levels(variable))
  }

  render_raw_table <- function(name){
    renderTable({
      dat <- get_raw_data()
      if(is.null(dat))
        return(NULL)
      v <- dat[[input[[name]]]]
      table(v, useNA = "always", dnn=name)
    })
  }

  render_table <- function(name){
    renderTable({
      v <- eval(parse(text=paste0("get_", name)))()
      table(v, useNA = "always", dnn=name)
    })
  }

  get_weights <- reactive({
    get_numeric("weights")
  })

  get_hiv <-  reactive({
    get_logical("hiv")
  })

  get_undiagnosed <-  reactive({
    get_logical("undiagnosed")
  })

  get_recent <-  reactive({
    get_logical("recent")
  })

  get_ever_test <-  reactive({
    get_logical("ever_test")
  })

  get_low_viral <-  reactive({
    get_logical("low_viral")
  })

  get_last_test <-  reactive({
    get_numeric("last_test")
  })

  get_treated <-  reactive({
    get_logical("treated")
  })

  get_rep_weights <- reactive({
    dat <- get_raw_data()
    rw_names <- input$rep_weights
    if(is.null(dat) || length(rw_names) == 0)
      return(NULL)
    rw <- lapply(as.list(rw_names), function(name){
      as.numeric(dat[[name]])
    })
    as.data.frame(rw)
  })

  get_last_test_upper <-  reactive({
    low <- get_numeric("last_test_lower")
    up <- get_numeric("last_test_upper")
    up[is.na(up) & !is.na(low) & low == max(low, na.rm=TRUE)]
    up
  })

  get_strata <-  reactive({
    get_categorical("strata")
  })


  output$contents <- renderTable({
    if(is.null(get_raw_data()))
      return(NULL)

    head(get_raw_data())
  })

  output$table <- renderDataTable({
    if(is.null(get_raw_data()))
      return(NULL)

    head(get_raw_data())
  })



  output$last_test_plot <- renderPlot({
    lt <- get_last_test()
    if(is.null(lt))
      return(NULL)
    print(qplot(lt, bins = 30) + xlab("Last HIV Test"))
  })

  output$last_test_errors <- renderText({
    lt <- get_last_test()
    if(is.null(lt))
      return("")
    txt <- ""
    if(sum(is.na(lt)) > length(lt) / 2)
      txt <- paste(txt, "Error: More than half of values are missing",sep="\n")
    lt <- na.omit(lt)
    if(any(lt <= 0))
      txt <- paste(txt, "Error: Last HIV test values <= 0 detected",sep="\n")
    if(any(lt > 35 * 365))
      txt <- paste(txt, "Error: Last HIV test values > 35 years ago",sep="\n")
    txt
  })

  output$hiv_desc <- render_table("hiv")
  output$hiv_desc_raw <- render_raw_table("hiv")

  output$undiagnosed_desc <- render_table("undiagnosed")
  output$undiagnosed_desc_raw <- render_raw_table("undiagnosed")

  output$recent_desc <- render_table("recent")
  output$recent_desc_raw <- render_raw_table("recent")

  output$low_viral_desc <- render_table("low_viral")
  output$low_viral_desc_raw <- render_raw_table("low_viral")

  output$ever_test_desc <- render_table("ever_test")
  output$ever_test_desc_raw <- render_raw_table("ever_test")

  output$strata_desc <- render_table("strata")
  output$strata_desc_raw <- render_raw_table("strata")

  output$treated_desc <- render_table("treated")
  output$treated_desc_raw <- render_raw_table("treated")

  incidence <- reactiveVal()
  nclicks <- reactiveVal(0)
  output$inc_results <- renderTable({
    nclicks(0)

    # Rerun on UI change
    undiagnosed <- get_undiagnosed()
    recent <- get_recent()
    low_viral <- get_low_viral()
    hiv <- get_hiv()
    last_test <- get_last_test()
    ever_test <- get_ever_test()
    treated <- get_treated()
    if(input$rita_type != "Treatment + Viral Load (RITA2)")
      treated <- NULL

    strata <- get_strata()

    weights <- get_weights()
    if(is.null(weights))
      weights <- rep(1, length(hiv))

    if(is.null(get_raw_data()))
      stop("No Data Uploaded")
    if(is.null(undiagnosed))
      stop("Required Variable Not Specified: Undiagnosed")
    if(is.null(hiv))
      stop("Required Variable Not Specified: HIV Status")
    if(is.null(recent))
      stop("Required Variable Not Specified: Assay Recent")
    if(is.null(last_test))
      stop("Required Variable Not Specified: Last Test")
    if(is.null(low_viral))
      stop("Required Variable Not Specified: Low Viral Load")
    #sm_strata <- min(table(strata))
    #if(!is.null(strata) && sm_strata < 5){
    #  stop("Stratifying Variable Has Strata With Too Few Observations")
    #}
    tau <- input$tau
    frr <- input$frr
    test_history_population <- input$test_history_population
    median_ttt <- input$median_ttt * 365
    if(is.na(median_ttt) && !is.null(treated)){
      stop("Please specify a median time to treatment")
    }
    lambda <- log(2) / median_ttt
    treat_surv <- 1 - pexp(1:ceiling(365*tau), lambda)
    inc <- list()
    if(is.null(strata)){
      inc[["all"]] <- rita_incidence(recent=recent,
                                        undiagnosed = undiagnosed,
                                        ever_hiv_test = ever_test,
                                        low_viral = low_viral,
                                        hiv = hiv,
                                        tslt=last_test,
                                        weights=weights,
                                        tau=tau,
                                        frr=frr,
                                        test_history_population=test_history_population,
                                     treated = treated,
                                     treat_surv = treat_surv
                      )
      result <- inc[[1]]

    }else{
      lvls <- levels(strata)
      for(lv in lvls){
        inc[[lv]] <- try(rita_incidence(recent=recent[strata == lv],
                                    undiagnosed = undiagnosed[strata == lv],
                                    ever_hiv_test = ever_test[strata == lv],
                                    low_viral = low_viral[strata == lv],
                                    hiv = hiv[strata == lv],
                                    tslt=last_test[strata == lv],
                                    weights=weights[strata == lv],
                                    tau=tau,
                                    frr=frr,
                                    test_history_population=test_history_population,
                                    treated = treated[strata == lv],
                                    treat_surv = treat_surv
                        ))
      }
      result <- inc
      for(i in 1:length(inc)){
        result[[i]] <- cbind(names(result)[[i]], result[[i]])
        row.names(result[[i]]) <- NULL
      }
      names(result) <- NULL
      result <- do.call(rbind,result)
      colnames(result)[1] <- "strata"
    }
    incidence(inc)
    boot_result(NULL)
    as.data.frame(result)
  }, rownames=FALSE,
  width="400px", digits=4)

  interruptor <- AsyncInterruptor$new()




  boot_result <- reactiveVal()
  observeEvent(input$run,{
    if(nclicks() != 0){
      print("Already running")
      return(NULL)
    }
    nclicks(nclicks() + 1)
    boot_result(data.frame(Status="Running..."))
    if(is.null(incidence()))
      return(NULL)
    #incc <- incidence()
    rep_weights <- get_rep_weights()
    nrep <- ncol(rep_weights)


    # Rerun on UI change
    undiagnosed <- get_undiagnosed()
    recent <- get_recent()
    low_viral <- get_low_viral()
    hiv <- get_hiv()
    last_test <- get_last_test()
    ever_test <- get_ever_test()
    treated <- get_treated()
    if(input$rita_type != "Treatment + Viral Load (RITA2)")
      treated <- NULL

    strata <- get_strata()

    weights <- get_weights()
    if(is.null(weights))
      weights <- rep(1, length(hiv))

    if(is.null(get_raw_data()))
      stop("No Data Uploaded")
    if(is.null(undiagnosed))
      stop("Required Variable Not Specified: Undiagnosed")
    if(is.null(hiv))
      stop("Required Variable Not Specified: HIV Status")
    if(is.null(recent))
      stop("Required Variable Not Specified: Assay Recent")
    if(is.null(last_test))
      stop("Required Variable Not Specified: Last Test")
    if(is.null(low_viral))
      stop("Required Variable Not Specified: Low Viral Load")
    #sm_strata <- min(table(strata))
    #if(!is.null(strata) && sm_strata < 5){
    #  stop("Stratifying Variable Has Strata With Too Few Observations")
    #}
    tau <- input$tau
    frr <- input$frr
    test_history_population <- input$test_history_population
    type <- input$type
    median_ttt <- input$median_ttt * 365
    if(is.na(median_ttt) && !is.null(treated)){
      stop("Please specify a median time to treatment")
    }
    lambda <- log(2) / median_ttt
    treat_surv <- 1 - pexp(1:ceiling(365*tau), lambda)
    #progress <- AsyncProgress$new(message="Generating Boostrap Samples")
    prog <- function(i, strata, nstrata){
      interruptor$execInterrupts()
      #progress$set(( (strata - 1) * nrep + i) / (nrep*nstrata))
    }
    inc <- list()
    if(is.null(strata)){
      callback <- function(k) prog(k, 1, 1)
      inc[["all"]] <- rita_bootstrap(recent=recent,
                                     undiagnosed = undiagnosed,
                                     ever_hiv_test = ever_test,
                                     low_viral = low_viral,
                                     hiv = hiv,
                                     tslt=last_test,
                                     weights=weights,
                                     tau=tau,
                                     frr=frr,
                                     rep_weights = rep_weights,
                                     rep_weight_type = type,
                                     test_history_population=test_history_population,
                                     treated = treated,
                                     treat_surv = treat_surv,
                                     show_progress = callback
      )
      inc[[1]] <- cbind(
        quantity = row.names(inc[[1]]),
        inc[[1]]
      )
      result <- inc[[1]]

    }else{
      lvls <- levels(strata)
      nstrata <- length(lvls)
      for(i in 1:length(lvls)){
        lv <- lvls[i]
        callback <- function(k) prog(k, i, nstrata)
        inc[[lv]] <- try(rita_bootstrap(recent=recent[strata == lv],
                                        undiagnosed = undiagnosed[strata == lv],
                                        ever_hiv_test = ever_test[strata == lv],
                                        low_viral = low_viral[strata == lv],
                                        hiv = hiv[strata == lv],
                                        tslt=last_test[strata == lv],
                                        weights=weights[strata == lv],
                                        tau=tau,
                                        frr=frr,
                                        rep_weights = rep_weights[strata == lv,],
                                        rep_weight_type = type,
                                        test_history_population=test_history_population,
                                        treated = treated[strata == lv],
                                        treat_surv = treat_surv,
                                        show_progress = callback
        ))
        inc[[lv]] <- cbind(
          quantity = row.names(inc[[lv]]),
          stratum = lv,
          inc[[lv]]
        )
      }
      result <- dplyr::bind_rows(inc) %>% dplyr::arrange(quantity,stratum)
      row.names(result) <- NULL
    }
    boot_result(as.data.frame(result))

    # progress <- AsyncProgress$new(message="Generating Boostrap Samples")
    # prog <- function(i, strata, nstrata){
    #   interruptor$execInterrupts()
    #   progress$set(( (strata - 1) * nrep + i) / (nrep*nstrata))
    # }
    # stratified <- length(incc) > 1
    # result <- finally(
    #   catch(
    #     future({
    #       blist <- list()
    #       nstrata <- length(incc)
    #       for(i in 1:nstrata){
    #         callback <- function(k) prog(k, i, nstrata)
    #         if(is.null(rep_weights))
    #           boot <- as.data.frame(summary(bootstrap_incidence(incc[[i]],
    #                                                             nrep=nrep,
    #                                                             show_progress=callback)))
    #         else
    #           boot <- as.data.frame(summary(bootstrap_incidence(incc[[i]],
    #                                                             rep_weights=rep_weights,
    #                                                             type=type,
    #                                                             show_progress=callback)))
    #         if(stratified){
    #           boot <- cbind(names(incc)[i], boot)
    #           names(boot)[1] <- "strata"
    #         }
    #         if(nrow(boot) > 1){
    #           boot <- cbind(row.names(boot), boot)
    #           names(boot)[1] <- "age_subgroup"
    #         }
    #         blist[[i]] <- boot
    #       }
    #       do.call(rbind, blist)
    #     })  %...>% boot_result,
    #     function(e) {
    #       boot_result(NULL)
    #       print(e$message)
    #       showNotification(e$message)
    #     }
    #   ),
    #   function(){
    #     progress$sequentialClose()
    #     nclicks(0)
    #   }
    # )

    NULL
  })
  output$bootstrap <- renderTable({
    shiny::req(boot_result())
  },digits=4)

  observeEvent(input$cancel,{
    print("cancel")
    interruptor$interrupt("User Interrupt")
  })

  output$download_example <- downloadHandler(
    filename = "assay_data.csv",
    content = function(file) {
      #loc <- paste0(system.file("shiny_ui", package = "rita"),"/assay_data.csv")
      loc <- "/usr/local/lib/R/site-library/rita/shiny_ui/assay_data.csv"
      file.copy(loc, file)
    }
  )

})
