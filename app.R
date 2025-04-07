# Load config
source("src/config.R")

# Load plotting functions
source("src/functions/plots.R")
source("src/functions/shiny_plots.R")

# Load file loc config
config <- yaml::yaml.load_file('wbuild.yaml')
ANNOTATION <- config$ANNOTATION
RAW_HPO <- config$RAW_HPO
RAW_RNA <- config$RAW_RNA
RAW_Protein <- config$RAW_Protein
PROC_DATA <- config$PROC_DATA

# Load sample annotation
sa <- fread(ANNOTATION)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]
sa[ , solved_class := "unsolved" ]
sa[  CATEGORY %in% c("IIb", "IIc" ), solved_class := "candidate" ]
sa[ CATEGORY %in% c("I.m", "I", "IIa", "III" ), solved_class := "solved" ]

sa2 <- fread("../data/omicsDiagnostics/proteomics_annotation.tsv")
sa2 <- sa2[SAMPLE_ID %in% sa$SAMPLE_ID  ]
sa2 <- sa2[USE_FOR_PROTEOMICS_PAPER == T]
sa2 <- sa2[!is.na(exome_ID)]


SSS <- readRDS("/Users/Mitya/Desktop/2024_latest_SCRIPTS_data/inegr/hpo_sss/SemSim_MWES_PP.rds") %>% as.data.table() 
SSS <- SSS[ type == "gene"]
SSSx <- SSS[PHENOME_ID %in% sa2$exome_ID ]
uniqueN(SSSx$PHENOME_ID)
colnames(SSSx)
SSSx <- SSSx[ , c("PHENOME_ID", "geneID"  , "SemSim" )]
setnames(SSSx, "PHENOME_ID", "exome_ID")
SSSx <- SSSx[!is.na(geneID)]

exomiser <- fread("/Users/Mitya/Desktop/2024_latest_SCRIPTS_data/inegr/exomiser/PP_exomiser.csv")
uniqueN(exomiser$exome_ID)
colnames(exomiser)
exomiser <- exomiser[exome_ID %in% sa2$exome_ID  , c("exome_ID" , "geneID", "EXOMISER_GENE_VARIANT_SCORE" ,"EXOMISER_GENE_COMBINED_SCORE", "Pval" , "EXOMISER_GENE_PHENO_SCORE")]
exomiser <- exomiser[!is.na(geneID)]
exomiser[ , geneID := toupper(geneID)]
uniqueN(exomiser$exome_ID)

exomiserX <- merge(exomiser, SSSx, by = c( "exome_ID" , "geneID"), all = T  )
exomiserX[ is.na(SemSim), SemSim := EXOMISER_GENE_PHENO_SCORE]
exomiserX[is.na(exomiserX)] <- 0
sa2 <- sa2[exome_ID %in% unique(exomiserX$exome_ID) ]

exomiserXo <- merge(sa2[, c("SAMPLE_ID", "exome_ID") ], exomiserX, by = "exome_ID")
uniqueN(exomiserXo$SAMPLE_ID)

exomiserXo <- exomiserXo[!duplicated(exomiserXo)]

# Read integrated omics file
rp <- readRDS(paste0(PROC_DATA, "/integration/patient_omics_full.RDS")) %>% as.data.table()
rp[ is.na(Semantic_sim), Semantic_sim:= 0]

#rp <- rp[outlier_class != 'non_outlier' | causal_gene == T ] #subset significant
#rp <- rp[ !is.na(PROTEIN_ZSCORE ) & !is.na(RNA_ZSCORE )  ]
rp <- rp[SAMPLE_ID %in% unique(exomiserXo$SAMPLE_ID)]
uniqueN(rp$SAMPLE_ID)
rpX <- merge(rp, exomiserXo, by = c( "SAMPLE_ID" , "geneID"), all.x = T  )
rp <- rpX
colnames(rp)
rp[ is.na(SemSim), SemSim:= 0]
rp$EXOMISER_GENE_PHENO_SCORE <- NULL
rp[ is.na(EXOMISER_GENE_VARIANT_SCORE), EXOMISER_GENE_VARIANT_SCORE:= 0]
rp[ is.na(EXOMISER_GENE_COMBINED_SCORE), EXOMISER_GENE_COMBINED_SCORE:= 0]
rp[ is.na(Pval), Pval:= 1]
rp[, sample_gene := paste0(SAMPLE_ID, "_", geneID)]

# Subset amples in the annotation
sa <- sa[SAMPLE_ID %in% unique(rp$SAMPLE_ID)]
 

complexes <- readRDS(paste0(PROC_DATA, "/Complexes/Complex_outliers_LIMMA.rds")) %>% as.data.table()
 
# Subset diagnosed cases and candidates
complexes <- complexes[SAMPLE_ID %in% sa$SAMPLE_ID]



# Load patient's HPO 
pat_hpo <- fread(RAW_HPO)

# Subset diagnosed cases and candidates
pat_hpo <- pat_hpo[SAMPLE_ID %in% sa$SAMPLE_ID, c("SAMPLE_ID", "HPO_term")]
pat_hpo <-pat_hpo[!duplicated(pat_hpo ),]

# Load patient's affected organ systems
patients <- fread(paste0(PROC_DATA, '/HPO/Patients_affected_organs.tsv'))
# Subset diagnosed cases and candidates
patients <- patients[SAMPLE_ID %in% sa$SAMPLE_ID]





 


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("omicsDiagnostics APP"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel( width = 2,
                      
            # Solved class check box 
            checkboxGroupInput("solvedF", 
                              "Diagnosis filter", 
                               choices = list("Unsolved" = "unsolved" , 
                                                        "Candidate" = "candidate", 
                                                        "Diagnosed" = "solved"),
                               selected = c( "unsolved", "candidate", "solved"  )), # "no rare variant" ,                       
            
            # Select sample           
            selectInput(inputId = "Sample",
                        label = "Choose Sample ID:", # 
                        multiple = T,
                        choices = NULL),
  
            
            
            # Horizontal line ----
            tags$hr(),
            
            # RNA Z slider
            sliderInput(inputId = "z_RNA", label ="|RNA z-score|:", min = 0, round(max(abs(rp[!is.na(RNA_ZSCORE)]$RNA_ZSCORE)),0) , 0),
            
            # Protein Z slider
            sliderInput(inputId = "z_Protein", label ="|Protein z-score|:",  min =  0, round(max(abs(rp[!is.na(PROTEIN_ZSCORE)]$PROTEIN_ZSCORE)),0) , 2),
            
            # Horizontal line ----
            tags$hr(),
            
            # SemSim slider
            sliderInput(inputId = "SemSim", label = "Phenotype similarity %:",  min =  0, round(max(abs(rp$SemSim)), 0), 1),
            
            # Exomiser slider
            # Checkbox what to plot 
            
            # Horizontal line ----
            tags$hr(),
            
            
            # Variant checkbox 
            checkboxGroupInput("variantsF", 
                               "Variant filter", 
                               choices = list("No variant" = "no rare variant" , 
                                              "1 RV" = "1 rare variant", 
                                              "2+ RV" = "rare pot. biallelic variants"),
                               selected = c( "1 rare variant", "rare pot. biallelic variants"  )), # "no rare variant" , 
            
            # Horizontal line ----
            tags$hr(),
            
            # Gene selector 
            textInput("genesF", "Select genes of interest", "ACAD9"), # A character string of the text input.
            helpText("Use comma spase ', ' \nto specify more than one"),
            
            # selectInput(inputId = "genesF",
            #             label = "Choose Sample ID:", # 
            #             multiple = T,
            #             choices = unique(rp[order(geneID)]$geneID)),
            
            # Horizontal line ----
            tags$hr(),
            h6("Test version Dima")
        ),

        # Show a plot of the generated distribution
        mainPanel(width = 10,
          h1("RNA and Protein integration"),
          tableOutput("patientInfo"),
          
          # first row 
          fluidRow(
            
            
            #HPO info and fig
            column(4,
                   tabsetPanel(type = "tabs",
                               tabPanel("Anatogram", plotOutput("anatogram", height = "600px") ) ,
                               tabPanel("Phenotypes",  tableOutput("patientHPO" ) ) 
                               ), # end tabset 
                   ), #End column HPO
            
            
            # OMICs plots
            column(8,
          tabsetPanel(type = "tabs",
                      tabPanel("Integration", plotlyOutput("patientPlot", height = "600px" ) ), # ,  width = "600px", height = "600px" 
                      tabPanel("Protein volcano", plotlyOutput("proteinVolcanoPlot", height = "600px") ),  # , width = "600px", height = "600px"  
                      tabPanel("RNA volcano" , plotlyOutput("rnaVolcanoPlot", height = "600px" ) )  # , width = "600px", height = "600px" 
                      ), # end tabset for plots
            ), #End column OMICS
          
          
          ), # end fluid row 
          # Horizontal line ----
          tags$hr(),
          
          # Second row 
          fluidRow(
            
            # Protein complexes
            column(4,
                   plotlyOutput("complexVolcanoPlot") # , width = "500px", height = "500px"
            ), #End column Protein complexes
            
            
            # Rank RNA 
            column(4,
                   plotlyOutput("rankRNAPlot" ) # , width = "500px", height = "500px"
            ), #End column Rank RNA 
            
            # Rank Protein 
            column(4,
                   plotlyOutput("rankProteinPlot" ) # , width = "500px", height = "500px"
            ), #End column Rank Protein 
            
            
      
            
       
            
          ), # end fluid row 2 
          
          # Horizontal line ----
          tags$hr(),
          
          
          h2("OMICS table"),
          # last row 
          DTOutput("omicsTable",  width = "100%")

        )
    )
)


server <- function(input,output, session){
  
  # Reactive expression for filtered samples based on solved class
  filtered_samples <- reactive({
    sa[solved_class %in% input$solvedF, SAMPLE_ID]
  })
  
  # Update the sample selection input based on the filtered samples
  observe({
    updateSelectInput(session, "Sample", choices = filtered_samples())
  })
  

  output$patientInfo <- renderTable({ 
    sa[(SAMPLE_ID %in% input$Sample)  , 
       c("SAMPLE_ID", "gender", "KNOWN_MUTATION", "CATEGORY", "solved_class", "TISSUE", "PROTEOMICS_BATCH") ]
  })
  
  output$patientHPO <- renderTable({ 
    pat_hpo[(SAMPLE_ID %in% input$Sample)  , 
       c("HPO_term") ]
  })
  
  # Reactive data filtering
  pat <- reactive({
    tmp <- rp[abs(RNA_ZSCORE) >= input$z_RNA &
                abs(PROTEIN_ZSCORE) >= input$z_Protein &
                SemSim >= input$SemSim &
                gene_class %in% input$variantsF]
    tmp <- tmp[SAMPLE_ID %in% sa[solved_class %in% input$solvedF]$SAMPLE_ID]
    # if(input$genesF != ""){
    #   tmp <- tmp[geneID %in% toupper(unlist(strsplit(input$genesF, ", ")))]
    # }
    if(!is.null(input$Sample)){
      tmp <- tmp[(SAMPLE_ID %in% input$Sample)]
    }
    tmp
  })
  
  output$patientPlot <- renderPlotly({
    patientOM <- pat()
    plot_patient_integration_plotly(patientOM)
  })
  
  output$proteinVolcanoPlot <- renderPlotly({
    data <- pat()
    plot <- plot_protein_volcano_plotly(data)
    ggplotly(plot)
  })
  
  output$rnaVolcanoPlot <- renderPlotly({
    data <- pat()
    plot <- plot_rna_volcano_plotly(data)
    ggplotly(plot)
  })
  
  output$anatogram <- renderPlot({
    req(input$Sample)
    plot_anatogram(input$Sample[1], patients, sa)
  })
  
  output$complexVolcanoPlot <- renderPlotly({
    req(input$Sample)
    plot_complex_volcano_plotly(input$Sample[1], complexes)
  })
  
  
  output$rankRNAPlot <- renderPlotly({
    req(input$genesF)
    gene <- toupper(unlist(strsplit(input$genesF, ", ")))[1]
    plot_rank_rna_plotly(rp, gene)
  })
  
  output$rankProteinPlot <- renderPlotly({
    req(input$genesF)
    gene <- toupper(unlist(strsplit(input$genesF, ", ")))[1]
    plot_rank_protein_plotly(rp, gene)
  })
  
  # Render variant table
  output$omicsTable <- DT::renderDT({ 
    selected_columns <- c("SAMPLE_ID", "geneID","rare", "potential_biallelic",
                          "EXOMISER_GENE_COMBINED_SCORE","SemSim",  "Pval",  "PROTEIN_ZSCORE",
                          "EXOMISER_GENE_VARIANT_SCORE", "MOI", "MOI_match","OMIM" ,
                          "PROTEIN_FC", "RNA_ZSCORE", "RNA_FC",
                          "outlier_class", "validated",
                            "RNA_PVALUE", "RNA_PADJ", 
                           "PROTEIN_PVALUE", "PROTEIN_PADJ"  )
    pat()[, ..selected_columns]
  })


}

colnames(rp)

# Run the application 
shinyApp(ui = ui, server = server)



