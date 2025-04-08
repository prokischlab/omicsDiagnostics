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

# Load configuration with error handling
tryCatch({
  config <- yaml::yaml.load_file('omicsDiagnosticsAPP/config.yaml')
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

# ... rest of the app.R file remains unchanged ... 