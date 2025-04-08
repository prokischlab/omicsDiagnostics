# Install required packages if not already installed
if (!require("rsconnect")) install.packages("rsconnect")

# Define required packages
required_packages <- c(
  # Core Shiny packages
  "shiny", "shinyjs", "shinybusy", "shinyWidgets", "shinythemes", "bslib",
  # Data handling
  "data.table", "yaml",
  # Visualization
  "plotly", "ggplot2", "gganatogram", "DT",
  # UI enhancements
  "tippy"
)

# Install required packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Authorize your account
rsconnect::setAccountInfo(
  name = "prokischlab",
  token = "C02A55A3009C9FDF18E35C24FCFE7F40",
  secret = "<SECRET>"
)

# Function to deploy with retry logic
deploy_with_retry <- function(max_attempts = 3, wait_time = 60) {
  for (attempt in 1:max_attempts) {
    tryCatch({
      message(sprintf("Deployment attempt %d of %d", attempt, max_attempts))
      
      # First, verify all required files exist
      required_files <- c(
        "app.R",
        "src/functions/shiny_plots.R",
        "requirements.txt",
        "omicsDiagnosticsAPP/sample_annotation.RDS",
        "omicsDiagnosticsAPP/patient_omics.RDS",
        "omicsDiagnosticsAPP/patient_hpo.RDS",
        "omicsDiagnosticsAPP/patients_organs.RDS",
        "omicsDiagnosticsAPP/complexes.RDS",
        "omicsDiagnosticsAPP/config.yaml"
      )
      
      missing_files <- required_files[!file.exists(required_files)]
      if (length(missing_files) > 0) {
        stop("Missing required files: ", paste(missing_files, collapse = ", "))
      }
      
      rsconnect::deployApp(
        appDir = ".",
        appName = "omicsDiagnostics",
        appTitle = "omicsDiagnostics App",
        account = "prokischlab",
        forceUpdate = TRUE,
        appFiles = required_files,
        logLevel = "verbose"
      )
      
      # If we get here, deployment was successful
      message("Deployment successful!")
      return(TRUE)
      
    }, error = function(e) {
      message(sprintf("Attempt %d failed: %s", attempt, e$message))
      
      if (attempt < max_attempts) {
        message(sprintf("Waiting %d seconds before retrying...", wait_time))
        Sys.sleep(wait_time)
      } else {
        message("All deployment attempts failed. Please try again later or contact support.")
        return(FALSE)
      }
    })
  }
}

# Run deployment with retry logic
deploy_with_retry() 