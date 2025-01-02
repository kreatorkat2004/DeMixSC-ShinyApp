devtools::install_github('wwylab/DeMixSC')

library(DeMixSC)
library(plotly)
library(dplyr)
library(shinyjs)

server <- function(input, output, session) {
  # Reactive values
  values <- reactiveValues(
    data_uploaded = FALSE,
    analysis_complete = FALSE,
    results = NULL
  )
  
  # File upload handlers
  observeEvent(input$countsFile, {
    req(input$countsFile)
    values$counts <- read.csv(input$countsFile$datapath, row.names=1)
    values$data_uploaded <- !is.null(values$counts) && !is.null(values$annotations)
  })
  
  observeEvent(input$annotFile, {
    req(input$annotFile)
    values$annotations <- read.csv(input$annotFile$datapath)
    values$data_uploaded <- !is.null(values$counts) && !is.null(values$annotations)
  })
  
  # Custom Analysis
  observeEvent(input$runAnalysis, {
    req(values$counts)
    withProgress(message = 'Running DeMixSC analysis...', value = 0, {
      incProgress(0.1, detail = "Loading data...")
      
      try({
        values$results <- DeMixSC(
          option = "user.defined",
          benchmark.mode = FALSE,
          mat.target = values$counts,
          min.expression = input$minExpression,
          scale.factor = input$scaleFactor
        )
        
        values$analysis_complete <- TRUE
        incProgress(1)
      })
    })
  })
  
  # Package Data Analysis - Retina
  observeEvent(input$runRetina, {
    withProgress(message = 'Loading Retina data...', value = 0, {
      if(input$retinaDataset == "AMD Cohort") {
        amd_cohort_raw_counts = get("retina_amd_cohort_raw_counts", envir = asNamespace("DeMixSC"))
        amd_cohort_clinical_info = get("retina_amd_cohort_clinical_info", envir = asNamespace("DeMixSC"))
        
        values$results <- DeMixSC(
          option = "retina",
          benchmark.mode = FALSE,
          mat.target = amd_cohort_raw_counts,
          min.expression = 3,
          scale.factor = 1e5
        )
      }
      values$analysis_complete <- TRUE
    })
  })
  
  # Package Data Analysis - HGSC
  observeEvent(input$runHGSC, {
    withProgress(message = 'Loading HGSC data...', value = 0, {
      if(input$hgscDataset == "Lee Cohort") {
        lee_cohort_raw_counts = get("hgsc_lee_cohort_raw_counts", envir = asNamespace("DeMixSC"))
        lee_cohort_clinical_info = get("hgsc_lee_cohort_clinical_info", envir = asNamespace("DeMixSC"))
        
        values$results <- DeMixSC(
          option = "hgsc",
          benchmark.mode = FALSE,
          mat.target = lee_cohort_raw_counts,
          min.expression = 5,
          scale.factor = 1e6
        )
      }
      values$analysis_complete <- TRUE
    })
  })
  
  # Plot outputs
  output$plot1 <- renderPlotly({
    req(values$results)
    plot_data <- values$results$cell.type.proportions
    
    plot_ly(type = "bar", x = rownames(plot_data), y = as.numeric(plot_data[,1])) %>%
      layout(title = "Cell Type Proportions",
             xaxis = list(title = "Cell Types"),
             yaxis = list(title = "Proportion"))
  })
  
  output$caption1 <- renderUI({
    req(values$results)
    HTML("<p>Cell type proportions estimated from deconvolution analysis.</p>")
  })
  
  output$plot2 <- renderPlotly({
    req(values$results)
    # Add second visualization based on results
  })
  
  output$caption2 <- renderUI({
    req(values$results)
    HTML("<p>Additional visualization of deconvolution results.</p>")
  })
  
  # Results download handler
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste("demixsc-results-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(values$results$cell.type.proportions, file)
    }
  )
  
  # Plot output area
  output$plotOutputArea <- renderUI({
    req(values$analysis_complete)
    div(
      h4("Results"),
      DTOutput("resultsTable")
    )
  })
  
  output$resultsTable <- renderDT({
    req(values$results)
    datatable(values$results$cell.type.proportions,
              options = list(scrollX = TRUE))
  })
}
