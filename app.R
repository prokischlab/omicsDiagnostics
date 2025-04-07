# Load required libraries
library(shiny)
library(data.table)
library(plotly)
library(DT)
library(yaml)
library(gganatogram)
library(shinyjs)
library(shinybusy)

# Load config and functions
source("src/config.R")
source("src/functions/plots.R")
source("src/functions/shiny_plots.R")

# Load configuration with error handling
tryCatch({
  config <- yaml::yaml.load_file('wbuild.yaml')
  ANNOTATION <- config$ANNOTATION
  RAW_HPO <- config$RAW_HPO
  RAW_RNA <- config$RAW_RNA
  RAW_Protein <- config$RAW_Protein
  PROC_DATA <- config$PROC_DATA
  
  # Validate paths
  if (!all(file.exists(c(ANNOTATION, RAW_HPO, RAW_RNA, RAW_Protein)))) {
    stop("One or more required data files are missing")
  }
}, error = function(e) {
  stop("Error loading configuration: ", e$message)
})

# Load sample annotation with error handling
tryCatch({
  sa <- fread(ANNOTATION)
  sa <- sa[USE_FOR_PROTEOMICS_PAPER == TRUE]
  sa[, solved_class := "unsolved"]
  sa[CATEGORY %in% c("IIb", "IIc"), solved_class := "candidate"]
  sa[CATEGORY %in% c("I.m", "I", "IIa", "III"), solved_class := "solved"]
}, error = function(e) {
  stop("Error loading sample annotation: ", e$message)
})

# Load additional data with error handling
tryCatch({
  sa2 <- fread("../data/omicsDiagnostics/proteomics_annotation.tsv")
  sa2 <- sa2[SAMPLE_ID %in% sa$SAMPLE_ID]
  sa2 <- sa2[USE_FOR_PROTEOMICS_PAPER == TRUE]
  sa2 <- sa2[!is.na(exome_ID)]
  
  SSS <- readRDS("/Users/Mitya/Desktop/2024_latest_SCRIPTS_data/inegr/hpo_sss/SemSim_MWES_PP.rds") %>% as.data.table()
  SSS <- SSS[type == "gene"]
  SSSx <- SSS[PHENOME_ID %in% sa2$exome_ID]
  SSSx <- SSSx[, c("PHENOME_ID", "geneID", "SemSim")]
  setnames(SSSx, "PHENOME_ID", "exome_ID")
  SSSx <- SSSx[!is.na(geneID)]
  
  exomiser <- fread("/Users/Mitya/Desktop/2024_latest_SCRIPTS_data/inegr/exomiser/PP_exomiser.csv")
  exomiser <- exomiser[exome_ID %in% sa2$exome_ID, 
                      c("exome_ID", "geneID", "EXOMISER_GENE_VARIANT_SCORE", 
                        "EXOMISER_GENE_COMBINED_SCORE", "Pval", "EXOMISER_GENE_PHENO_SCORE")]
  exomiser <- exomiser[!is.na(geneID)]
  exomiser[, geneID := toupper(geneID)]
  
  exomiserX <- merge(exomiser, SSSx, by = c("exome_ID", "geneID"), all = TRUE)
  exomiserX[is.na(SemSim), SemSim := EXOMISER_GENE_PHENO_SCORE]
  exomiserX[is.na(exomiserX)] <- 0
  sa2 <- sa2[exome_ID %in% unique(exomiserX$exome_ID)]
  
  exomiserXo <- merge(sa2[, c("SAMPLE_ID", "exome_ID")], exomiserX, by = "exome_ID")
  exomiserXo <- exomiserXo[!duplicated(exomiserXo)]
}, error = function(e) {
  stop("Error loading additional data: ", e$message)
})

# Load integrated omics data with error handling
tryCatch({
  rp <- readRDS(paste0(PROC_DATA, "/integration/patient_omics_full.RDS")) %>% as.data.table()
  rp[is.na(Semantic_sim), Semantic_sim := 0]
  rp <- rp[SAMPLE_ID %in% unique(exomiserXo$SAMPLE_ID)]
  rpX <- merge(rp, exomiserXo, by = c("SAMPLE_ID", "geneID"), all.x = TRUE)
  rp <- rpX
  
  # Clean and prepare data
  rp[is.na(SemSim), SemSim := 0]
  rp$EXOMISER_GENE_PHENO_SCORE <- NULL
  rp[is.na(EXOMISER_GENE_VARIANT_SCORE), EXOMISER_GENE_VARIANT_SCORE := 0]
  rp[is.na(EXOMISER_GENE_COMBINED_SCORE), EXOMISER_GENE_COMBINED_SCORE := 0]
  rp[is.na(Pval), Pval := 1]
  rp[, sample_gene := paste0(SAMPLE_ID, "_", geneID)]
  
  # Subset samples in the annotation
  sa <- sa[SAMPLE_ID %in% unique(rp$SAMPLE_ID)]
}, error = function(e) {
  stop("Error loading integrated omics data: ", e$message)
})

# Load additional data with error handling
tryCatch({
  complexes <- readRDS(paste0(PROC_DATA, "/Complexes/Complex_outliers_LIMMA.rds")) %>% as.data.table()
  complexes <- complexes[SAMPLE_ID %in% sa$SAMPLE_ID]
  # Load patient's HPO 
  pat_hpo <- fread(RAW_HPO)
  
  hpo <- get_ontology("datasets/hp.obo", extract_tags="everything")
  hpo_terms <- data.table(
    HPO_ID = hpo$id,
    HPO_term = hpo$name
  )

  
  pat_hpo <- merge(pat_hpo,hpo_terms, by = "HPO_ID" )
  
  
  # Subset diagnosed cases and candidates
  pat_hpo <- pat_hpo[SAMPLE_ID %in% sa$SAMPLE_ID, c("SAMPLE_ID", "HPO_term")]
  pat_hpo <-pat_hpo[!duplicated(pat_hpo ),]
  
  rm(hpo, hpo_terms)
  
  
  pat_hpo <- pat_hpo[SAMPLE_ID %in% sa$SAMPLE_ID, c("SAMPLE_ID", "HPO_term")]
  pat_hpo <- pat_hpo[!duplicated(pat_hpo)]
  
  patients <- fread(paste0(PROC_DATA, '/HPO/Patients_affected_organs.tsv'))
  patients <- patients[SAMPLE_ID %in% sa$SAMPLE_ID]
}, error = function(e) {
  stop("Error loading additional data: ", e$message)
})

# Define UI for application
ui <- fluidPage(
  titlePanel("omicsDiagnostics APP"),
  
  # Add loading spinner
  shinyjs::useShinyjs(),
  shinybusy::add_busy_spinner(spin = "fading-circle"),
  
  sidebarLayout(
    sidebarPanel(
      width = 2,
      
      # Diagnosis filter
      h4("Diagnosis Filter"),
      checkboxGroupInput("solvedF",
                        label = NULL,
                        choices = list("Unsolved" = "unsolved",
                                     "Candidate" = "candidate",
                                     "Diagnosed" = "solved"),
                        selected = c("unsolved", "candidate", "solved")),
      
      # Sample selection
      h4("Sample Selection"),
      selectInput(inputId = "Sample",
                 label = "Choose Sample ID:",
                 multiple = TRUE,
                 choices = NULL),
      
      hr(),
      
      # Z-score filters
      h4("Z-score Filters"),
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
      
      hr(),
      
      # Phenotype similarity
      h4("Phenotype Similarity"),
      sliderInput(inputId = "SemSim",
                 label = "Phenotype similarity %:",
                 min = 0,
                 max = round(max(abs(rp$SemSim)), 0),
                 value = 1),
      
      hr(),
      
      # Variant filter
      h4("Variant Filter"),
      checkboxGroupInput("variantsF",
                        label = NULL,
                        choices = list("No variant" = "no rare variant",
                                     "1 RV" = "1 rare variant",
                                     "2+ RV" = "rare pot. biallelic variants"),
                        selected = c("1 rare variant", "rare pot. biallelic variants")),
      
      hr(),
      
      # Gene selection
      h4("Gene Selection"),
      textInput("genesF",
                "Select genes of interest",
                "ACAD9"),
      helpText("Use comma space ', ' to specify multiple genes"),
      
      hr(),
      
      # Version info
      h6("Test version Dima")
    ),
    
    mainPanel(
      width = 10,
      h1("RNA and Protein Integration"),
      
      # Patient info table
      DTOutput("patientInfo"),
      
      # First row of plots
      fluidRow(
        column(4,
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Anatogram", plotOutput("anatogram", height = "600px")),
                 tabPanel("Phenotypes", DTOutput("patientHPO"))
               )
        ),
        column(8,
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Integration", plotlyOutput("patientPlot", height = "600px")),
                 tabPanel("Protein volcano", plotlyOutput("proteinVolcanoPlot", height = "600px")),
                 tabPanel("RNA volcano", plotlyOutput("rnaVolcanoPlot", height = "600px"))
               )
        )
      ),
      
      hr(),
      
      # Second row of plots
      fluidRow(
        column(4, plotlyOutput("complexVolcanoPlot")),
        column(4, plotlyOutput("rankRNAPlot")),
        column(4, plotlyOutput("rankProteinPlot"))
      ),
      
      hr(),
      
      # OMICS table
      h2("OMICS Table"),
      DTOutput("omicsTable", width = "100%")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Add error handling for reactive expressions
  filtered_samples <- reactive({
    tryCatch({
      sa[solved_class %in% input$solvedF, SAMPLE_ID]
    }, error = function(e) {
      showNotification("Error filtering samples", type = "error")
      return(NULL)
    })
  })
  
  # Update sample selection
  observe({
    tryCatch({
      updateSelectInput(session, "Sample", choices = filtered_samples())
    }, error = function(e) {
      showNotification("Error updating sample selection", type = "error")
    })
  })
  
  # Patient info table
  output$patientInfo <- renderDT({
    tryCatch({
      datatable(
        sa[SAMPLE_ID %in% input$Sample,
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
        pat_hpo[SAMPLE_ID %in% input$Sample, c("HPO_term")],
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
      tmp <- rp[abs(RNA_ZSCORE) >= input$z_RNA &
                 abs(PROTEIN_ZSCORE) >= input$z_Protein &
                 SemSim >= input$SemSim &
                 gene_class %in% input$variantsF]
      tmp <- tmp[SAMPLE_ID %in% sa[solved_class %in% input$solvedF]$SAMPLE_ID]
      if (!is.null(input$Sample)) {
        tmp <- tmp[SAMPLE_ID %in% input$Sample]
      }
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
      if (is.null(patientOM)) return(NULL)
      plot_patient_integration_plotly(patientOM)
    }, error = function(e) {
      showNotification("Error rendering integration plot", type = "error")
      return(NULL)
    })
  })
  
  output$proteinVolcanoPlot <- renderPlotly({
    tryCatch({
      data <- pat()
      if (is.null(data)) return(NULL)
      plot <- plot_protein_volcano_plotly(data)
      ggplotly(plot)
    }, error = function(e) {
      showNotification("Error rendering protein volcano plot", type = "error")
      return(NULL)
    })
  })
  
  output$rnaVolcanoPlot <- renderPlotly({
    tryCatch({
      data <- pat()
      if (is.null(data)) return(NULL)
      plot <- plot_rna_volcano_plotly(data)
      ggplotly(plot)
    }, error = function(e) {
      showNotification("Error rendering RNA volcano plot", type = "error")
      return(NULL)
    })
  })
  
  output$anatogram <- renderPlot({
    tryCatch({
      req(input$Sample)
      plot_anatogram(input$Sample[1], patients, sa)
    }, error = function(e) {
      showNotification("Error rendering anatogram", type = "error")
      return(NULL)
    })
  })
  
  output$complexVolcanoPlot <- renderPlotly({
    tryCatch({
      req(input$Sample)
      plot_complex_volcano_plotly(input$Sample[1], complexes)
    }, error = function(e) {
      showNotification("Error rendering complex volcano plot", type = "error")
      return(NULL)
    })
  })
  
  output$rankRNAPlot <- renderPlotly({
    tryCatch({
      req(input$genesF)
      gene <- toupper(unlist(strsplit(input$genesF, ", ")))[1]
      plot_rank_rna_plotly(rp, gene)
    }, error = function(e) {
      showNotification("Error rendering RNA rank plot", type = "error")
      return(NULL)
    })
  })
  
  output$rankProteinPlot <- renderPlotly({
    tryCatch({
      req(input$genesF)
      gene <- toupper(unlist(strsplit(input$genesF, ", ")))[1]
      plot_rank_protein_plotly(rp, gene)
    }, error = function(e) {
      showNotification("Error rendering protein rank plot", type = "error")
      return(NULL)
    })
  })
  
  # OMICS table with error handling
  output$omicsTable <- renderDT({
    tryCatch({
      selected_columns <- c("SAMPLE_ID", "geneID", "rare", "potential_biallelic",
                          "EXOMISER_GENE_COMBINED_SCORE", "SemSim", "Pval", "PROTEIN_ZSCORE",
                          "EXOMISER_GENE_VARIANT_SCORE", "MOI", "MOI_match", "OMIM",
                          "PROTEIN_FC", "RNA_ZSCORE", "RNA_FC",
                          "outlier_class", "validated",
                          "RNA_PVALUE", "RNA_PADJ",
                          "PROTEIN_PVALUE", "PROTEIN_PADJ")
      
      datatable(
        pat()[, ..selected_columns],
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
        ),
        extensions = 'Buttons'
      )
    }, error = function(e) {
      showNotification("Error rendering OMICS table", type = "error")
      return(NULL)
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)



