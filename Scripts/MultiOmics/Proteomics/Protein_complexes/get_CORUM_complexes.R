#'---
#' title: Get protein complexes from CORUM 
#' author: smirnovd
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  output:
#'  - corum: '`sm config["PROC_DATA"] + "/Complexes/CORUM.tsv"`'
#'  type: script
#'---

source(snakemake@input$config)
#Download CORUM all complexes
corum <- fread("curl http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip | funzip")

corum <- corum[Organism == 'Human'] # Subset for human data
corum <- corum[, c("ComplexName",  "subunits(Gene name)")]
colnames(corum) <- c("COMPLEX"  , "geneID")

corum[COMPLEX %like% "Ecsit complex" , geneID := paste0(geneID , ";ACAD9")] # Add known member of Ecsit complex

# Format data
corum <- corum[!duplicated(corum ) ]
corum <- as.data.table(separate_rows(corum, geneID))
corum <- corum[!duplicated(corum ) ]
corum <- corum[! (geneID %in% c('', ' ', "1", "A" ) )]
corum[, N_subunits:= .N, by = .(COMPLEX)]

# Use complexes with at least two subunits
corum <-  corum[N_subunits >= 2 ]
corum <- corum[, c("geneID", "COMPLEX", "N_subunits"  )]
corum <-  corum[!duplicated(corum ) ]
corum$geneID <-  toupper(corum$geneID)

paste("Number of unique protein complexes:", uniqueN(corum$COMPLEX) )
paste("median number of subunits:", median( corum[ , .N, by = COMPLEX]$N) )
paste("maximum number of subunits:", max( corum[ , .N, by = COMPLEX]$N) )

write_tsv(corum,  snakemake@output$corum)
# write_tsv(corum, '/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/Complexes/CORUM.tsv')




