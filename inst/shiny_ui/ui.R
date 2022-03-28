library(shiny)
library(shinyWidgets)
# Define UI for application that draws a histogram
shinyUI(
  navbarPage('RITA Incidence',
             #tag$script(HTML("updateFirstMultiInput = function(x){ var event = new Event('input'); $('.multi-wrapper .search-input').get(0).dispatchEvent(event);};Shiny.addCustomMessageHandler('updateFirstMultiInput', updateFirstMultiInput);")),
             tabPanel('Introduction',
                      #includeScript("updateMulti.js"),
                      h4('Welcome to The RITA Incidence Estimator'),
                      br(),
                      p("The purpose of this tool is the estimation of incidence from a cross-sectional survey utilizing a recency assay."),
                      br(),
                      h5('Please proceed to the', em('Load Data'), 'tab')
             ),
             tabPanel('Data',
                      fluidPage(
                        # Application title
                        titlePanel("Load Data"),

                        # Sidebar with a slider input for number of bins
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("file1", "Choose CSV File",
                                      accept = c(
                                        "text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")
                            ),
                            conditionalPanel("output.table != null",
                                             selectizeInput("hiv", "HIV+:",
                                                            c("")),
                                             selectizeInput("undiagnosed", "Undiagnosed:",
                                                            c("")),
                                             selectizeInput("ever_test", "Ever Had An HIV Test:",
                                                            c("")),
                                             selectizeInput("low_viral", "Low/Undetectable Viral Load:",
                                                            c("")),
                                             selectizeInput("recent", "Assay Recent:",
                                                            c("")),
                                             selectizeInput("last_test", "Time Since Last HIV Test (Years):",
                                                            c("")),
                                             selectizeInput("strata", "Stratify By (Optional):",
                                                            c("")),
                                             tags$hr(),
                                             selectizeInput("weights", "Weights (Optional):",
                                                            c("")),
                                             pickerInput(inputId = "rep_weights",
                                                         label = "Replication Weights (Optional):",
                                                         choices = c(""),
                                                         multiple=TRUE,
                                                         options = list(`actions-box` = TRUE,
                                                                        `live-search`=TRUE,
                                                                        `none-selected-text`="Choose Variable")
                                             )
                            ),
                            width=5),
                          mainPanel(
                            conditionalPanel("input.hiv != \"\"",
                                             h2("HIV Descriptives:"),
                                             p("Raw Values:"),
                                             tableOutput("hiv_desc_raw"),
                                             p("Processed:"),
                                             tableOutput("hiv_desc")
                            ),
                            conditionalPanel("input.undiagnosed != \"\"",
                                             h2("Undiagnosed Descriptives:"),
                                             p("Raw Values:"),
                                             tableOutput("undiagnosed_desc_raw"),
                                             p("Processed:"),
                                             tableOutput("undiagnosed_desc")
                            ),
                            conditionalPanel("input.ever_test != \"\"",
                                             h2("Undiagnosed Descriptives:"),
                                             p("Raw Values:"),
                                             tableOutput("ever_test_desc_raw"),
                                             p("Processed:"),
                                             tableOutput("ever_test_desc")
                            ),
                            conditionalPanel("input.low_viral != \"\"",
                                             h2("Low Viral Load Descriptives:"),
                                             p("Raw Values:"),
                                             tableOutput("low_viral_desc_raw"),
                                             p("Processed:"),
                                             tableOutput("low_viral_desc")
                            ),
                            conditionalPanel("input.recent != \"\"",
                                             h2("Assay Recent Descriptives:"),
                                             p("Raw Values:"),
                                             tableOutput("recent_desc_raw"),
                                             p("Processed:"),
                                             tableOutput("recent_desc")
                            ),
                            conditionalPanel("input.strata != \"\"",
                                             h2("Strata Descriptives:"),
                                             p("Raw Values:"),
                                             tableOutput("strata_raw"),
                                             p("Processed:"),
                                             tableOutput("strata_desc")
                            ),
                            conditionalPanel("input.last_test != \"\"",
                                             h2("Last Test Descriptives:"),
                                             plotOutput("last_test_plot", width = "100%", height = "200px"),
                                             verbatimTextOutput("last_test_errors")
                            ),
                            dataTableOutput("table"),
                            width = 7
                          )
                        )
                      )
             ),
             tabPanel('Analysis',
                      fluidPage(
                        # Application title
                        titlePanel("Analysis"),
                        # Sidebar with a slider input for number of bins
                        sidebarLayout(
                          sidebarPanel(
                              numericInput("tau","Recency Period (tau)", min = 0, value = 2),
                              numericInput("frr","Reference False Recency Rate (FRR)", min = 0, value =  .0055),
                              selectizeInput("test_history_population", "Testing History Population", choices = c("undiagnosed","negative"), selected="undiagnosed"),
                              h3("Bootstrap Confidence Intervals"),
                              selectInput("type",
                                          "Replicate Weight Type:",
                                          choices = c("Jackknife"="JK1","Bootstrap"="bootstrap","BRR", "Fay"),
                                          selected="JK1"
                              ),
                              actionButton('run', 'Run'),
                              actionButton('cancel', 'Cancel')
                          ),
                          mainPanel(
                            h3("Incidence Results"),
                            tableOutput("inc_results"),
                            #h3("Bootstrap Intervals"),
                            conditionalPanel("(input.rep_weights != null) || (input.design_clusters != \"\") || (input.design_strata != \"\")",
                                             h3("Survey Bootstrap Intervals")
                            ),
                            conditionalPanel("(input.rep_weights == null) && (input.design_clusters == \"\") && (input.design_strata == \"\")",
                                             h3("Bootstrap Intervals")
                            ),
                            tableOutput("bootstrap")
                          )
                        )
                      )
             )
  )
)
