suppressPackageStartupMessages({
  stringsAsFactors=F
# Load required libraries
library(shiny)
library(data.table)
library(plotly)
library(DT)
library(yaml)
library(gganatogram)
library(shinyjs)
library(shinybusy)
library(shinyWidgets)
library(shinythemes)
library(tippy)
library(bslib)
})

# Define color and shape scales
outlier_colors <- c(
  "non_outlier" = "gray80",
  "RNA-only" = "blue",
  "protein-only" = "red",
  "RNA-and-protein" = "purple"
)

gene_shapes <- c(
  "no data" = 20,
  "candidate" = 17,
  "solved" = 15
)

# Load plotting functions
source("shiny_plots.R")

# Outlier threshold definition 
PADJ_THRESHOLD =  0.1
# Outlier threshold definition for rna-and-protein outliers
ZSCORE_THRESHOLD = 3
LOG2FC_THRESHOLD = 1

MIN_SAMPLE_NUM_PROT <- 0.52 #0.76 # % of samples protein not detected in


# Define colours for plotting 

variant_colors <- c("non_coding" = "gray81",  
                    "splice" = "darkolivegreen2" ,  
                    "stop" =  "aquamarine3",  
                    "coding" = "salmon1",   
                    "frameshift" = "lightskyblue3",  
                    "synonymous" = "moccasin",
                    "no rare variant" = "white",
                    "rare" = "white", 
                    "rare " = "white", 
                    "potential_biallelic" = "gray70",
                    "potential biallelic SemSim > 2" = "gray50",
                    "potential biallelic + Semantic sim > 2" = "gray50")

variant_colors2 <- c("no rare variant" = "white",
                     "SemSim > 2"= "gray90",
                     "rare" = "gray70", 
                     "rare SemSim > 2" = "gray50",
                     "potential_biallelic" = "gray30",
                     "potential biallelic SemSim > 2" = "gray10") 


gene_size <- c("no data" = 2.5, 
               "no rare variant" = 2.5, 
               "1 rare variant" = 3.1,
               "rare pot. biallelic variants" = 4.5) 

text_color <- c("no data" = "grey20", 
                "no rare variant" = "grey20", 
                "1 rare variant" = "black", 
                "rare pot. biallelic variants" = "black")




# Load config and functions
# source("src/config.R")
# source("src/functions/plots.R")

# Load configuration with error handling
tryCatch({
  config <- yaml::yaml.load_file('config.yaml')
  ANNOTATION <- config$ANNOTATION
  PROC_DATA <- config$PROC_DATA
  
  # Validate paths
  if (!all(file.exists(c(ANNOTATION)))) {
    stop("Required data files are missing")
  }
}, error = function(e) {
  stop("Error loading configuration: ", e$message)
})

# Load sample annotation with error handling
tryCatch({
  sa <- readRDS(ANNOTATION)
}, error = function(e) {
  stop("Error loading sample annotation: ", e$message)
})

# Load integrated omics data with error handling
tryCatch({
  rp <- readRDS(paste0(PROC_DATA, "/patient_omics.RDS"))
  rp[ , SemSim := Semantic_sim_scaled]
}, error = function(e) {
  stop("Error loading integrated omics data: ", e$message)
})

# Load additional data with error handling
tryCatch({
  complexes <- readRDS(paste0(PROC_DATA, "/complexes.RDS"))
  pat_hpo <- readRDS(paste0(PROC_DATA, "/patient_hpo.RDS"))
  patients <- readRDS(paste0(PROC_DATA, "/patients_organs.RDS"))
}, error = function(e) {
  stop("Error loading additional data: ", e$message)
})

# Add data caching at the start
data_cache <- reactiveValues(
  rp = NULL,
  sa = NULL,
  pat_hpo = NULL,
  complexes = NULL,
  patients = NULL
)

# Load and cache data
observe({
  tryCatch({
    # Cache integrated omics data
    data_cache$rp <- rp
    
    # Cache sample annotation
    data_cache$sa <- sa
    
    # Cache HPO data
    data_cache$pat_hpo <- pat_hpo
    
    # Cache complexes data
    data_cache$complexes <- complexes
    
    # Cache patients data
    data_cache$patients <- patients
  }, error = function(e) {
    showNotification("Error caching data", type = "error")
  })
})

# Define UI for application
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly"),
  
  titlePanel("omicsDiagnostics APP"),
  
  # Add loading spinner and progress bar
  shinyjs::useShinyjs(),
  shinybusy::add_busy_spinner(spin = "fading-circle"),
  shinybusy::add_busy_bar(color = "#FF0000"),
  
  sidebarLayout(
    sidebarPanel(
      width = 2,
      
      # Diagnosis and Sample Selection
      wellPanel(
        h4("Diagnosis Filter"),
        checkboxGroupInput("solvedF",
                          label = NULL,
                          choices = list("Candidate" = "candidate",
                                       "Diagnosed" = "solved"),
                          selected = c("candidate", "solved")),
        
        h4("Sample Selection"),
        pickerInput(inputId = "Sample",
                    label = "Choose Sample ID:",
                    multiple = TRUE,
                    choices = NULL,
                    options = list(
                      `live-search` = TRUE,
                      `actions-box` = TRUE
                    ))
      ),
      
      # Filters
      wellPanel(
        h4("Filters"),
        
        # Z-score filters
        sliderInput(inputId = "z_RNA",
                    label = "|RNA z-score|:",
                    min = 0,
                    max = round(max(abs(rp[!is.na(RNA_ZSCORE)]$RNA_ZSCORE)), 0),
                    value = 0),
        
        sliderInput(inputId = "z_Protein",
                    label = "|Protein z-score|:",
                    min = 0,
                    max = round(max(abs(rp[!is.na(PROTEIN_ZSCORE)]$PROTEIN_ZSCORE)), 0),
                    value = 2),
        
        # Phenotype similarity
        sliderInput(inputId = "SemSim",
                    label = "Phenotype similarity %:",
                    min = 0,
                    max = 100,
                    value = 1),
        
        # Variant filter
        checkboxGroupInput("variantsF",
                          label = "Variant Filter",
                          choices = list("No variant" = "no rare variant",
                                       "1 RV" = "1 rare variant",
                                       "2+ RV" = "rare pot. biallelic variants"),
                          selected = c("1 rare variant", "rare pot. biallelic variants")),
        
        # Gene Selection
        textInput("genesF",
                  "Select genes of interest",
                  "ACAD9"),
        helpText("Use comma space ', ' to specify multiple genes")
      ),
      
      # Version info
      h6("Developed by Dima Smirnov")
    ),
    
    mainPanel(
      width = 10,
      h1("RNA and Protein Integration"),
      
      # Add summary statistics
      fluidRow(
        column(12,
               wellPanel(
                 h4("Summary Statistics"),
                 textOutput("summaryStats")
               )
        )
      ),
      
      # Patient info table with download button
      fluidRow(
        column(12,
               h3("Patient Information"),
               downloadButton("downloadPatientInfo", "Download Patient Info"),
               DTOutput("patientInfo")
        )
      ),
      
      # First row of plots with download buttons
      fluidRow(
        column(4,
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Anatogram", 
                         downloadButton("downloadAnatogram", "Download Anatogram"),
                         plotOutput("anatogram", height = "600px")),
                 tabPanel("Phenotypes", 
                         downloadButton("downloadHPO", "Download HPO Data"),
                         DTOutput("patientHPO"))
               )
        ),
        column(8,
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Integration", 
                         downloadButton("downloadIntegration", "Download Integration Plot"),
                         plotlyOutput("patientPlot", height = "600px")),
                 tabPanel("Protein volcano", 
                         downloadButton("downloadProteinVolcano", "Download Protein Volcano"),
                         plotlyOutput("proteinVolcanoPlot", height = "600px")),
                 tabPanel("RNA volcano", 
                         downloadButton("downloadRNAVolcano", "Download RNA Volcano"),
                         plotlyOutput("rnaVolcanoPlot", height = "600px"))
               )
        )
      ),
      
      hr(),
      
      # Second row of plots with download buttons
      fluidRow(
        column(4, 
               downloadButton("downloadComplexVolcano", "Download Complex Volcano"),
               plotlyOutput("complexVolcanoPlot")),
        column(4, 
               downloadButton("downloadRankRNA", "Download RNA Rank"),
               plotlyOutput("rankRNAPlot")),
        column(4, 
               downloadButton("downloadRankProtein", "Download Protein Rank"),
               plotlyOutput("rankProteinPlot"))
      ),
      
      hr(),
      
      # OMICS table with enhanced features
      h2("OMICS Table"),
      downloadButton("downloadOMICSTable", "Download OMICS Table"),
      DTOutput("omicsTable", width = "100%")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Cache filtered data
  filtered_data <- reactive({
    tryCatch({
      req(data_cache$rp)
      tmp <- data_cache$rp[abs(RNA_ZSCORE) >= input$z_RNA &
                           abs(PROTEIN_ZSCORE) >= input$z_Protein &
                           Semantic_sim_scaled >= input$SemSim &
                           gene_class %in% input$variantsF]
      tmp <- tmp[SAMPLE_ID %in% data_cache$sa[solved_class %in% input$solvedF]$SAMPLE_ID]
      if (!is.null(input$Sample)) {
        tmp <- tmp[SAMPLE_ID %in% input$Sample]
      }
      tmp
    }, error = function(e) {
      showNotification("Error filtering data", type = "error")
      return(NULL)
    })
  })
  
  # Update sample selection
  observe({
    req(data_cache$sa)
    updatePickerInput(session, "Sample", 
                     choices = data_cache$sa[solved_class %in% input$solvedF, SAMPLE_ID])
  })
  
  # Summary statistics
  output$summaryStats <- renderText({
    data <- filtered_data()
    if (is.null(data)) return("No data available")
    
    paste(
      "Total samples:", length(unique(data$SAMPLE_ID)), "\n",
      "Total genes:", length(unique(data$geneID)), "\n",
      "Outliers:", sum(data$outlier_class != "non_outlier"), "\n",
      "RNA outliers:", sum(data$outlier_class == "RNA-only"), "\n",
      "Protein outliers:", sum(data$outlier_class == "protein-only"), "\n",
      "RNA and Protein outliers:", sum(data$outlier_class == "RNA-and-protein")
    )
  })
  
  # Patient info table
  output$patientInfo <- renderDT({
    tryCatch({
      datatable(
        data_cache$sa[SAMPLE_ID %in% input$Sample,
           c("SAMPLE_ID", "gender", "KNOWN_MUTATION", "CATEGORY", "solved_class", "TISSUE", "PROTEOMICS_BATCH")],
        options = list(dom = 't', pageLength = 5)
      )
    }, error = function(e) {
      showNotification("Error rendering patient info", type = "error")
      return(NULL)
    })
  })
  
  # Patient HPO table
  output$patientHPO <- renderDT({
    tryCatch({
      datatable(
        data_cache$pat_hpo[SAMPLE_ID %in% input$Sample, c("HPO_term")],
        options = list(dom = 't', pageLength = 5)
      )
    }, error = function(e) {
      showNotification("Error rendering HPO data", type = "error")
      return(NULL)
    })
  })
  
  # Filtered data
  pat <- reactive({
    tryCatch({
      tmp <- filtered_data()
      if (is.null(tmp)) return(NULL)
      tmp
    }, error = function(e) {
      showNotification("Error filtering data", type = "error")
      return(NULL)
    })
  })
  
  # Plot outputs with error handling
  output$patientPlot <- renderPlotly({
    tryCatch({
      patientOM <- pat()
      if (is.null(patientOM)) {
        showNotification("No data available for plotting", type = "warning")
        return(NULL)
      }
      if (nrow(patientOM) == 0) {
        showNotification("No data points available after filtering", type = "warning")
        return(NULL)
      }
      plot_patient_integration_plotly(patientOM)
    }, error = function(e) {
      showNotification(paste("Error rendering integration plot:", e$message), type = "error")
      return(NULL)
    })
  })
  
  output$proteinVolcanoPlot <- renderPlotly({
    tryCatch({
      data <- pat()
      if (is.null(data)) {
        showNotification("No data available for plotting", type = "warning")
        return(NULL)
      }
      if (nrow(data) == 0) {
        showNotification("No data points available after filtering", type = "warning")
        return(NULL)
      }
      plot_protein_volcano_plotly(data)
    }, error = function(e) {
      showNotification(paste("Error rendering protein volcano plot:", e$message), type = "error")
      return(NULL)
    })
  })
  
  output$rnaVolcanoPlot <- renderPlotly({
    tryCatch({
      data <- pat()
      if (is.null(data)) {
        showNotification("No data available for plotting", type = "warning")
        return(NULL)
      }
      if (nrow(data) == 0) {
        showNotification("No data points available after filtering", type = "warning")
        return(NULL)
      }
      plot_rna_volcano_plotly(data)
    }, error = function(e) {
      showNotification(paste("Error rendering RNA volcano plot:", e$message), type = "error")
      return(NULL)
    })
  })
  
  output$anatogram <- renderPlot({
    tryCatch({
      req(input$Sample)
      plot_anatogram(input$Sample[1], data_cache$patients, data_cache$sa)
    }, error = function(e) {
      showNotification("Error rendering anatogram", type = "error")
      return(NULL)
    })
  })
  
  output$complexVolcanoPlot <- renderPlotly({
    tryCatch({
      req(input$Sample)
      plot_complex_volcano_plotly(input$Sample[1], data_cache$complexes)
    }, error = function(e) {
      showNotification("Error rendering complex volcano plot", type = "error")
      return(NULL)
    })
  })
  
  output$rankRNAPlot <- renderPlotly({
    tryCatch({
      req(input$genesF)
      gene <- toupper(unlist(strsplit(input$genesF, ", ")))[1]
      plot_rank_rna_plotly(data_cache$rp, gene)
    }, error = function(e) {
      showNotification("Error rendering RNA rank plot", type = "error")
      return(NULL)
    })
  })
  
  output$rankProteinPlot <- renderPlotly({
    tryCatch({
      req(input$genesF)
      gene <- toupper(unlist(strsplit(input$genesF, ", ")))[1]
      plot_rank_protein_plotly(data_cache$rp, gene)
    }, error = function(e) {
      showNotification("Error rendering protein rank plot", type = "error")
      return(NULL)
    })
  })
  
  # OMICS table with error handling
  output$omicsTable <- renderDT({
    tryCatch({
      data <- filtered_data()
      if (is.null(data)) return(NULL)
      
      selected_columns <- c("SAMPLE_ID", "geneID", "rare", "potential_biallelic",
                           "Semantic_sim_scaled", "PROTEIN_ZSCORE", "RNA_ZSCORE",
                           "MOI", "MOI_match", "OMIM",
                           "PROTEIN_FC", "RNA_FC",
                           "outlier_class", "validated",
                           "RNA_PVALUE", "RNA_PADJ",
                           "PROTEIN_PVALUE", "PROTEIN_PADJ")
      
      # Format numeric columns for better readability
      data[, Semantic_sim_scaled := round(Semantic_sim_scaled, 2)]
      data[, PROTEIN_ZSCORE := round(PROTEIN_ZSCORE, 2)]
      data[, RNA_ZSCORE := round(RNA_ZSCORE, 2)]
      data[, PROTEIN_FC := round(PROTEIN_FC, 2)]
      data[, RNA_FC := round(RNA_FC, 2)]
      data[, RNA_PVALUE := format(RNA_PVALUE, scientific = TRUE, digits = 2)]
      data[, RNA_PADJ := format(RNA_PADJ, scientific = TRUE, digits = 2)]
      data[, PROTEIN_PVALUE := format(PROTEIN_PVALUE, scientific = TRUE, digits = 2)]
      data[, PROTEIN_PADJ := format(PROTEIN_PADJ, scientific = TRUE, digits = 2)]
      
      # Convert logical columns to Yes/No for better filtering
      data[, rare := ifelse(rare, "Yes", "No")]
      data[, potential_biallelic := ifelse(potential_biallelic, "Yes", "No")]
      data[, validated := ifelse(validated, "Yes", "No")]
      data[, MOI_match := ifelse(MOI_match, "Yes", "No")]
      
      datatable(
        data[, ..selected_columns],
        options = list(
          scrollX = TRUE,
          pageLength = 25,
          lengthMenu = c(10, 25, 50, 100),
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
          columnDefs = list(
            list(targets = c(2, 3, 12, 13), # Columns with Yes/No values
                 searchable = TRUE,
                 searchPanes = list(show = TRUE)),
            list(targets = c(4, 5, 6, 10, 11), # Numeric columns
                 searchable = TRUE,
                 type = "num"),
            list(targets = c(0, 1), # Sample and gene IDs
                 searchable = TRUE,
                 searchPanes = list(show = TRUE)),
            list(targets = c(7, 8, 9), # MOI and OMIM
                 searchable = TRUE,
                 searchPanes = list(show = TRUE))
          ),
          searchPanes = list(
            layout = 'columns-4',
            threshold = 0.1
          ),
          order = list(list(4, 'desc')) # Default sort by semantic similarity
        ),
        extensions = c('Buttons', 'SearchPanes'),
        selection = 'none',
        filter = 'top',
        rownames = FALSE,
        class = 'display cell-border stripe',
        server = FALSE  # Use client-side processing
      ) %>%
        formatStyle(
          'outlier_class',
          target = 'row',
          backgroundColor = styleEqual(
            c("non_outlier", "RNA-only", "protein-only", "RNA-and-protein"),
            c('#FFFFFF', '#FFE4E1', '#E0FFFF', '#F0FFF0')
          )
        )
    }, error = function(e) {
      showNotification("Error rendering OMICS table", type = "error")
      return(NULL)
    })
  })
  
  # Add download handlers
  output$downloadPatientInfo <- downloadHandler(
    filename = function() {
      paste("patient_info_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data_cache$sa[SAMPLE_ID %in% input$Sample, 
                  c("SAMPLE_ID", "gender", "KNOWN_MUTATION", "CATEGORY", 
                    "solved_class", "TISSUE", "PROTEOMICS_BATCH")], 
                file, row.names = FALSE)
    }
  )
  
  # Enhanced error handling with user feedback
  observeEvent(input$Sample, {
    if (length(input$Sample) == 0) {
      showNotification("Please select at least one sample", type = "warning")
    }
  })
  
  observeEvent(input$genesF, {
    if (input$genesF == "") {
      showNotification("Please enter at least one gene", type = "warning")
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)



