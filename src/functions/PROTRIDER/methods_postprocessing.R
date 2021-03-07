############################################
### post-processing stuff
############################################

### add mito property to gene_list
library(plyr)
library(tidyr)




### TODO refractor for limma column from other_scripts/add_mito_genes.R




add_hgnc_symbol = function(prot_sample_table, hgnc_file="/s/project/mitoMultiOmics/processed_proteomics/Dima/hgnc_mapping.tsv") {
  hgnc_mapping = fread(dima_file)
  hgnc_mapping = hgnc_mapping[,-c("gene_id","To","gene_name", "gene_type")]
  colnames(hgnc_mapping)[1] = 'protein_id'
  joined_se_pval = join(prot_sample_table, hgnc_mapping)
  return(joined_se_pval)
}

# 
# source("/data/nasif12/home_if12/loipfins/gitlab/genetic_diagnosis/Scripts/_functions/gene_annotation/add_gene_info_cols.R")
# add_disease_cols = function(prot_sample_table){
#   dt_pval = add_hans_class(prot_sample_table, gene_name_col = 'hgnc_symbol', return_all_info = FALSE)
#   dt_pval = add_omim_cols(dt_pval, gene_name_col = 'hgnc_symbol', return_all_info = F, pmim_link =F)
#   data.frame(dt_pval)
# }



### check if sample gene has rare variant
# check_mutation_gene = function(sample_prot_table) {
#   extra_col = sapply(sample_prot_table$sample_prot, function(x){
#     tryCatch( {
#       ex = subset(sample_prot_table, sample_prot == x)
#       print(ex$sample_prot)
#       
#       ex_gene = ex$hgnc_symbol
#       ex_exome = ex$exome_id
#       vt = readRDS(paste0('/s/project/mitoMultiOmics/raw_data/helmholtz/',ex_exome,'/exomicout/paired-endout/processedData/vep_anno_',ex_exome,'_uniq_dt.Rds') )
#       
#       out_row = subset(vt, hgncid==ex_gene)[, c("mstype",'gt', 'MAX_AF')]
#       # out_row_short = subset(out_row, MAX_AF < 0.2 | is.na(MAX_AF))
#       print(out_row_short)
#       print('#####')
#     }, error = function(cond) {
#       print(paste0('ERROR: ',x))
#     })
#   })
#   
# }








