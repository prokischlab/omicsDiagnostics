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
      
      # Set the working directory to the app directory
      app_dir <- "omicsDiagnosticsAPP"
      setwd(app_dir)
      
      # List of required files for deployment
      required_files <- c(
        "app.R",
        "shiny_plots.R",
        "config.yaml",
        "shiny_data/sample_annotation.RDS",
        "shiny_data/patient_omics.RDS",
        "shiny_data/patient_hpo.RDS",
        "shiny_data/patients_organs.RDS",
        "shiny_data/complexes.RDS"
      )
      
      # Get the current directory
      current_dir <- getwd()
      message("Current directory: ", current_dir)
      
      # List all files in the current directory
      message("Files in current directory:")
      print(list.files())
      
      # List all files in shiny_data directory
      if (dir.exists("shiny_data")) {
        message("Files in shiny_data directory:")
        print(list.files("shiny_data"))
      } else {
        message("shiny_data directory does not exist")
      }
      
      # Verify all required files exist
      missing_files <- required_files[!file.exists(required_files)]
      if (length(missing_files) > 0) {
        stop("Missing required files: ", paste(missing_files, collapse = ", "))
      }
      
      # Deploy the application with memory optimization
      rsconnect::deployApp(
        appDir = ".",
        appName = "omicsDiagnostics",
        appTitle = "omicsDiagnostics App",
        account = "prokischlab",
        forceUpdate = TRUE,
        appFiles = required_files,
        logLevel = "verbose",
        # Add memory optimization settings
        appPrimaryDoc = "app.R",
        appMode = "shiny",
        appStore = "shinyapps.io",
        appVisibility = "public",
        # Set memory limits
        appInstanceSize = "large",  # Use large instance size
        appInstanceCount = 1,       # Single instance
        appInstanceMemory = 1024    # 1GB memory
      )
      
      # Reset working directory
      setwd("..")
      
      message("Deployment successful!")
      return(TRUE)
      
    }, error = function(e) {
      # Reset working directory in case of error
      setwd("..")
      
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
