devtools::install_github('wwylab/DeMixSC')

library(shiny)
library(DeMixSC)

server <- function(input, output, session) {
  # Reactive values for both workflows
  customData <- reactiveValues(
    counts = NULL,
    annotations = NULL,
    results = NULL
  )
  
  packageData <- reactiveValues(
    counts = NULL,
    annotations = NULL,
    results = NULL
  )
  
  # Custom Data Handler
  observeEvent(input$countsFile, {
    req(input$countsFile)
    customData$counts <- read.csv(input$countsFile$datapath, row.names=1)
  })
  
  observeEvent(input$annotFile, {
    req(input$annotFile)
    customData$annotations <- read.csv(input$annotFile$datapath)
  })
  
  # Package Data Handler
  observeEvent(input$datasetChoice, {
    if(input$datasetChoice == "retina") {
      packageData$counts = get("retina_amd_cohort_raw_counts", 
                               envir = asNamespace("DeMixSC"))
      packageData$annotations = get("retina_amd_cohort_clinical_info", 
                                    envir = asNamespace("DeMixSC"))
    } else if(input$datasetChoice == "hgsc") {
      packageData$counts = get("hgsc_lee_cohort_raw_counts", 
                               envir = asNamespace("DeMixSC"))
      packageData$annotations = get("hgsc_lee_cohort_clinical_info", 
                                    envir = asNamespace("DeMixSC"))
    }
  })
  
  # Custom Analysis
  observeEvent(input$runCustomAnalysis, {
    req(customData$counts)
    withProgress(message = 'Running analysis...', {
      customData$results <- DeMixSC(
        option = "user.defined",
        benchmark.mode = FALSE,
        mat.target = customData$counts,
        min.expression = input$minExpression,
        scale.factor = input$scaleFactor
      )
    })
  })
  
  # Package Analysis
  observeEvent(input$runPackageAnalysis, {
    req(packageData$counts)
    withProgress(message = 'Running analysis...', {
      packageData$results <- DeMixSC(
        option = input$datasetChoice,
        benchmark.mode = FALSE,
        mat.target = packageData$counts,
        min.expression = input$minExpression,
        scale.factor = input$scaleFactor
      )
    })
  })
  
  # Output Rendering
  output$customPlot <- renderPlotly({
    req(customData$results)
    # Add visualization code here
  })
  
  output$packagePlot <- renderPlotly({
    req(packageData$results)
    # Add visualization code here
  })
  
  output$customResults <- renderDT({
    req(customData$results)
    datatable(customData$results$cell.type.proportions)
  })
  
  output$packageResults <- renderDT({
    req(packageData$results)
    datatable(packageData$results$cell.type.proportions)
  })
}
