# Integration of proteomics with genomics and transcriptomics increases the diagnosis rate of Mendelian disorders

This project contains different scripts to automatize and visualize analysis performed for the "Integration of proteomics with genomics and transcriptomics increases the diagnosis rate of Mendelian disorders" [paper](https://www.medrxiv.org/content/10.1101/2021.03.09.21253187v1).

[Webserver](https://prokischlab.github.io/omicsDiagnostics/#readme.html), produced as one of the outputs of the pipeline. 

## Project structure

This project is setup as a [wBuild workflow](https://github.com/gagneurlab/wBuild). This is an automatic build tool for R reports based on [snakemake](https://snakemake.readthedocs.io/en/stable/).

* The `wbuild.yaml` is the main configuration file to setup up the workflow
* The `Scripts` folder contains scripts which will be rendered as HTML reports
* The `src` folder contains additional helper functions and scripts
* The `Output` folder will contain all files produced in the analysis pipeline
    * `Output/html` contains the final HTML report

## Data and prerequisites 

This project depends on the packages `wBuild` and `PROTRIDER`

The pipeline starts with the series of files available via [Zenodo](https://zenodo.org/record/4501904): `DOI:10.5281/zenodo.4501904`

* `raw_data`
  * `proteomics_annotation.tsv` - sample annotation
  * `proteomics_not_normalized.tsv` - Proteomics intensity matrix
  * `raw_counts.tsv` - RNA-seq count matrix
  * `Patient_HPO_phenotypes.tsv` - Phenotype data recorded using HPO terms for diagnosed cases.
  * `enrichment_proportions_variants.tsv` - Results of rare variant enrichment/proportion analysis calculated on the full dataset.
  * `patient_variant_hpo_data.tsv` - Gene annotation for all individuals. 
    Since the genetic data are not publicly shareable, we provide only gene-level information for outlier genes only.


* `datasets`
  * `disease_genes.tsv` - List of Mendelian disease genes aggregated from several studies.
  * `HGNC_mito_groups.tsv` - Subset of [HGNC gene groups](https://www.genenames.org/tools/search/#!/groups?query=mitochondrial) related to mitochondria.
  
     Downloaded automatically:
  * `gencode.v29lift37.annotation.gtf.gz` - Gene-level model based on the [GENCODE 29 transcript model](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.annotation.gtf.gz) 
  * `Table_S1_gene_info_at_protein_level.xlsx` - Supplementary Tble1 from [GTEx proteomics study Jiang et al, 2020, Cell](https://www.cell.com/cell/fulltext/S0092-8674(20)31078-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867420310783%3Fshowall%3Dtrue)
  Data is available at the [GTEx page](https://storage.googleapis.com/gtex_egtex/proteomics/Table_S1_gene_info_at_protein_level.xlsx)
  * `allComplexes.txt` - CORUM protein complexes, available at CORUM [web page](http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip)

The proteomic raw data and MaxQuant search files have been deposited to the [ProteomeXchange Consortium](http://proteomecentral.proteomexchange.org) via the [PRIDE partner repository](https://www.ebi.ac.uk/pride/archive/projects/PXD022803) and can be accessed using the dataset identifier `PXD022803`

## Repository setup

First download the repo and its dependencies:

```
# analysis code
git clone https://github.com/prokischlab/omicsDiagnostics
cd omicsDiagnostics
```

and install wbuild using pip by running.

```
pip install wBuild
wBuild init
```

Since `wBuild init` will reset the current `Snakefile`, ` readme.md`, and `wbuild.yaml` we have to revert them again with git.

```
git checkout Snakefile
git checkout wbuild.yaml
git checkout readme.md
```


Specify correct file and folder locations in the `wbuild.yaml`.



## Install dependencies
* R packages
  * Make sure that `data.table` is installed or install with `install.packages("data.table")`
  * `Rscript src/installRPackages.R src/requirementsR.txt`

* Create Conda environment
  * `conda env create --name omicsDiagnostics --file=environment.yml`



## Run the full pipeline
To run the full pipeline, execute the following commands with 10 cores in parallel:

1) `conda activate omicsDiagnostics`

2) `snakemake graph`

3) `snakemake publish -j 10`

