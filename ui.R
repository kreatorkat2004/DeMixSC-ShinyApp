library(shiny)
library(plotly)
library(DT)

ui <- fluidPage(
  titlePanel("DeMixSC Analysis Interface"),
  
  tabsetPanel(
    # Tab 1: Custom Data Analysis
    tabPanel("Custom Data Analysis",
             sidebarLayout(
               sidebarPanel(
                 fileInput("countsFile", "Upload Count Matrix:", 
                           accept = c(".csv", ".txt")),
                 fileInput("annotFile", "Upload Cell Type Annotations:", 
                           accept = c(".csv", ".txt")),
                 numericInput("numClusters", "Number of Clusters:", 
                              value = 2, min = 2),
                 selectInput("normMethod", "Normalization Method:",
                             choices = c("LogNormalize", "SCTransform")),
                 actionButton("runCustomAnalysis", "Run Analysis")
               ),
               mainPanel(
                 plotlyOutput("customPlot"),
                 DTOutput("customResults")
               )
             )
    ),
    
    # Tab 2: Package Data Analysis
    tabPanel("Package Data Analysis",
             sidebarLayout(
               sidebarPanel(
                 selectInput("datasetChoice", "Select Dataset:",
                             choices = c("Retina" = "retina", 
                                         "HGSC" = "hgsc")),
                 numericInput("minExpression", "Minimum Expression:",
                              value = 3),
                 numericInput("scaleFactor", "Scale Factor:",
                              value = 1e5),
                 actionButton("runPackageAnalysis", "Run Analysis")
               ),
               mainPanel(
                 plotlyOutput("packagePlot"),
                 DTOutput("packageResults")
               )
             )
    )
  )
)
