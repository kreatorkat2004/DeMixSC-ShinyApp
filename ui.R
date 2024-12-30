library(shiny)
library(plotly)
library(shinyjs)
library(DT)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("DeMixSC: Deconvolution of Mixed-type Single Cell RNA-seq Data"),
  
  navbarPage("",
             # About tab
             tabPanel("About",
                      h4("Description"),
                      p("DeMixSC is a method for deconvoluting mixed-type scRNA-seq data..."),
                      p("Contact: [Your contact info]")
             ),
             
             # Analysis tab
             tabPanel("Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          # File upload section
                          fileInput("countsFile", "Upload Count Matrix:", accept = c(".csv", ".txt")),
                          fileInput("annotFile", "Upload Cell Type Annotations:", accept = c(".csv", ".txt")),
                          
                          # Parameters
                          numericInput("numClusters", "Number of Clusters:", value = 2, min = 2),
                          selectInput("normMethod", "Normalization Method:",
                                      choices = c("LogNormalize", "SCTransform")),
                          
                          # Action buttons
                          actionButton("runAnalysis", "Run Analysis"),
                          downloadButton("downloadResults", "Download Results")
                        ),
                        
                        mainPanel(
                          # Results tabs
                          tabsetPanel(
                            tabPanel("Visualization",
                                     plotlyOutput("umap"),
                                     plotlyOutput("proportions")
                            ),
                            tabPanel("Statistics",
                                     DTOutput("statsTable")
                            ),
                            tabPanel("Gene Expression",
                                     plotlyOutput("geneExpr")
                            )
                          )
                        )
                      )
             ),
             
             # Examples tab  
             tabPanel("Examples",
                      h4("Example Datasets"),
                      DTOutput("exampleTable")
             )
  )
)