if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager")
}

BiocManager::install("sva")
BiocManager::install("preprocessCore")
install.packages("doParallel")
install.packages("pbapply")
install.packages("nnls")

if (!"devtools" %in% rownames(installed.packages())) {
  install.packages('devtools')
}

devtools::install_github('wwylab/DeMixSC')

library(Seurat)
library(DeMixSC)

server <- function(input, output, session) {
  # Reactive values
  values <- reactiveValues(
    seurat_obj = NULL,
    results = NULL
  )
  
  # Load data
  observeEvent(input$countsFile, {
    req(input$countsFile)
    counts <- read.csv(input$countsFile$datapath)
    annot <- read.csv(input$annotFile$datapath)
    
    # Create Seurat object
    values$seurat_obj <- CreateSeuratObject(counts = counts)
  })
  
  # Run analysis
  observeEvent(input$runAnalysis, {
    req(values$seurat_obj)
    
    withProgress(message = 'Running DeMixSC analysis...', {
      # Normalize data
      values$seurat_obj <- NormalizeData(values$seurat_obj)
      
      # Find variable features
      values$seurat_obj <- FindVariableFeatures(values$seurat_obj)
      
      # Run DeMixSC
      values$results <- DeMixSC(values$seurat_obj,
                                num.clusters = input$numClusters)
    })
  })
  
  # Generate plots
  output$umap <- renderPlotly({
    req(values$results)
    # UMAP visualization code
  })
  
  output$proportions <- renderPlotly({
    req(values$results)
    # Cell type proportions plot
  })
  
  output$statsTable <- renderDT({
    req(values$results)
    # Statistics table
  })
  
  output$geneExpr <- renderPlotly({
    req(values$results)
    # Gene expression visualization
  })
  
  # Download handler
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste("demixsc-results-", Sys.Date(), ".zip", sep="")
    },
    content = function(file) {
      # Package results into zip file
    }
  )
}