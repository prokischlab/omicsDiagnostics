options(repos=structure(c(CRAN="https://cloud.r-project.org")), warn = -1)
suppressPackageStartupMessages(library(data.table))

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
  BiocManager::install("remotes")
}

install.packages("data.table")

args <- commandArgs(trailingOnly=TRUE)
packages <- fread(args[1], fill = TRUE)
packages <- packages[!startsWith(package, "#")]
packages$version <- as.character(packages$version)
packages[is.na(version), version := ""]
installed <- as.data.table(installed.packages())

for (pckg_name in unique(packages$package)) {
  package_dt <- packages[package == pckg_name]
  pckg_name <- tail(unlist(strsplit(pckg_name, split = "/")), n = 1)
  version <- package_dt$version
  
  if (pckg_name %in% installed$Package &
      (version == "" || installed[Package == pckg_name, Version] == version)
  ) {
    #message(paste(pckg_name, "already installed"))
  } else {
    if (package_dt$bioconductor == TRUE) {
      INSTALL <- BiocManager::install
    } else {
      INSTALL <- install.packages
    }
    package <- package_dt$package
    message(paste("install", package))
    INSTALL(package)
    message(paste("installed", package))
  }
}

if ("gganatogram" %in% installed$Package){
  message("gganatogram already installed")
} else {
  devtools::install_github("jespermaag/gganatogram")
}
  

if (!requireNamespace("OUTRIDER", quietly = TRUE) || packageVersion("OUTRIDER") < "1.20.1") {
  devtools::install_github("gagneurlab/OUTRIDER@1.20.1", dependencies = TRUE)
} else {
  message("OUTRIDER is already at version 1.20.1 or higher.")
}

options(warn = 0)