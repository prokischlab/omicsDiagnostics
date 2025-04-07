
sa <- fread("/data/agprokisch/working/dima/multiOMICs_integration/raw_data/proteomics_annotation.tsv")
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

sac <- fread("/data/agprokisch/working/dima/CowRD/data/project_data/proteome_analysis/raw_data/proteomics_annotation.tsv")
sac <- sac[USE_FOR_PROTEOMICS_PAPER == T, 1:9]
sa <- merge(sa, sac, by = "SAMPLE_ID")
rm(sac)



sss <- fread("/data/agprokisch/working/dima/CowRD/data/project_data/auto_analysis/semantic_sim.csv")
sss <- sss[PHENOME_ID %in% sa$exome_ID | PHENOME_ID %in% sa$GENOME_ID]
uniqueN(sss$PHENOME_ID)


exomiser <- fread("/data/agprokisch/working/dima/CowRD/data/project_data/auto_analysis/exomiser_table.csv")
setnames(exomiser, "#GENE_SYMBOL", "geneID")
exomiser <- exomiser[WES_ID %in%  sa$exome_ID]
exomiser <- exomiser[order(EXOMISER_GENE_COMBINED_SCORE, decreasing = T)]
exomiser[ , wes_gene := paste0(WES_ID, "_", geneID)]
exomiser <- exomiser[!duplicated(wes_gene)]
uniqueN(exomiser$WES_ID)

sse <- merge(exomiser, sss, by.x = c("WES_ID", "geneID"), by.y = c("PHENOME_ID", "geneID"), all = T )
colnames(sse)
sse <- sse[ , c("WES_ID", "geneID", "MOI_exomiser" , "Semantic_sim", "SemSim", "Rank_SemSim", "EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_GENE_VARIANT_SCORE")]
sse <- sse[!duplicated(sse)]

############################

out_path <- "/data/agprokisch/working/dima/CowRD/data/raw_data/WES_WGS/RV_tables/samples/"
vt <- c()
for ( wes in unique(c(sa$exome_ID)) ){
  out_wes <- paste0(out_path, wes, "/variants_acmg.txt")
  if (file.exists(out_wes)) {
    vt_a <- fread(out_wes)
    vt <- rbind(vt, vt_a)
  }
}
rm(vt_a)
uniqueN(vt$WES_ID) # 320


vt[vt == ""] <- NA
vt <- vt[!duplicated(vt)]

vt[!is.na(geneID_eva), geneID_eva := geneID_intervar]
vt <- vt[!is.na(geneID_eva)]
vt <- vt[geneID_eva != "."]
vt$geneID_intervar <- NULL
vt <- vt[!duplicated(vt)]

colnames(vt)
vt <- vt[ , c("WES_ID","geneID_eva","aux","gt", "N_var", "gene_class",  "MOI", 
              "ACMG", "clinvar_eva", "pLI", "oe_lof", "oe_mis", "obs_hom_lof", "obs_het_lof",
              "eva_class", "eva_func", "CADDraw", "CADDphred", "SIFT", 
              "wes_count_inhouse", "AF_inhouse", "AF_gnomad", "AF_1000G", "MAX_AF", "InterVar")]
vt <- vt[!duplicated(vt)]
vt <- vt[!duplicated(vt[ , c("WES_ID","geneID_eva","aux","gt", "N_var", "gene_class",  "MOI","ACMG", "CADDphred" )] ) ,]
vt <- vt[!is.na(geneID_eva)]

hgnc_map <-  getCurrentHumanMap() %>% as.data.table()
hgnc_map[ , Symbol := toupper(Symbol)]
hgnc_map[ , Approved.Symbol := toupper(Approved.Symbol)]
hgnc_map <- hgnc_map[!is.na(Symbol)]

vta <- merge(vt, hgnc_map, by.x = "geneID_eva", by.y = "Symbol", all.x = T, allow.cartesian=TRUE)
vta[!is.na(Approved.Symbol), Approved.Symbol := geneID_eva]
vta <- vta[!duplicated(vta)]
setnames(vta, "Approved.Symbol", "geneID")
vta <- vta[!is.na(geneID)]
vta[ACMG %in% c("VUS", "LB", "B") &  clinvar_eva %in% c("Likely_pathogenic", "Pathogenic"), ACMG := "LP"]

vts <- vta[ACMG %in% c("LP", "P"), c("WES_ID", "geneID_eva") ]
vts <- vts[, .N, by = c("WES_ID", "geneID_eva")]
setnames(vts, "N", "N_pathogenic")

vta <- merge(vta, vts, c("WES_ID", "geneID_eva"), all.x = T)
colnames(vta)
vts <- vta[, c("WES_ID", "geneID", "geneID_eva", "MOI", "N_var", "N_pathogenic", "gene_class",  "pLI" ) ]
vts <- vts[!duplicated(vts)]

rm(vt, hgnc_map, sss, exomiser)
########################################


rp <- readRDS("/data/agprokisch/working/dima/multiOMICs_integration/processed_data2/integration/patient_omics_full.RDS") %>% as.data.table()
rp <- readRDS("/data/agprokisch/working/dima/CowRD/data/project_data/proteome_analysis/processed_data/outrider_protrider_10_11.RDS") %>% as.data.table()

vts[geneID %in% unique(sse$geneID) | geneID %in% unique(rp$geneID) , gene := "other" ]
vts[geneID_eva %in% unique(sse$geneID) | geneID_eva %in% unique(rp$geneID), gene := "eva" ]
vts[gene == "eva", geneID :=  geneID_eva]
vts$geneID_eva <- NULL
vts$gene <- NULL

vtsse <- merge(vts, sse, by = c("WES_ID", "geneID"), all = T)
colnames(vtsse)
vtsse <- vtsse[ , c("WES_ID", "geneID", "gene_class", "N_pathogenic", "MOI_exomiser", "SemSim", "EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_GENE_VARIANT_SCORE" )]
setnames(vtsse, c("EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_GENE_VARIANT_SCORE"), c("exomiser", "exomiser_var_score"))
vtsse <- vtsse[!duplicated(vtsse)]

colnames(sa)
sa <- sa[, c("SAMPLE_ID", "exome_ID", "CATEGORY")]

colnames(rp)
rp <- merge( rp, sa,  by = c("SAMPLE_ID"))
rp <- rp[ , -c("rare", "potential_biallelic")]
sa[, solved := F]
sa[CATEGORY %in% c( "I", "IIa", "III") , solved := T]
setnames(rp, "gene_class", "gene_class_old")

rp <- merge(vtsse, rp, by.x = c("WES_ID", "geneID"), by.y = c("exome_ID", "geneID"), all.y = T)
rp[ , solved := F]
rp[CATEGORY %in% c( "I", "IIa", "III") , solved := T]
# rp[ WES_ID %in% sa[ solved == T]$exome_ID ,  solved := T]
rp <- rp[!duplicated(rp)]

rps <- rp[ (abs(PROTEIN_ZSCORE) >= 2) | ( is.na(PROTEIN_ZSCORE) & abs(RNA_ZSCORE) > 2 )]
rps <- rps[solved == F]
colnames(rps)
rps[is.na(gene_class),  gene_class := gene_class_old]
rps <- rps[  !is.na(gene_class) | !is.na(N_pathogenic) | !is.na(SemSim) | !is.na(exomiser)  ]
rps <- rps[ exomiser > 0.7 ]
rps <- rps[order( exomiser, PROTEIN_ZSCORE, RNA_ZSCORE, SemSim, decreasing = T)]
rps[exomiser_var_score > 0.7 & gene_class == "no rare variant", gene_class := "1 rare variant"  ]
rps[gene_class ==  "rare pot. biallelic variants", gene_class :=  "pot. biallelic"  ]
rps$MOI_exomiser <- NULL
rps$gene_class_old <- NULL

rpc <- rp[ causal_gene == T]
rps <- rps[ !(WES_ID %in% rpc$WES_ID) ]

rp_nj <- rps[!(  WES_ID %like% "JAP" ) ]
rp_nj <- rp_nj[order(WES_ID, exomiser, PROTEIN_ZSCORE, RNA_ZSCORE, SemSim, decreasing = T)]
rp_nj[ , .N, by = WES_ID]
uniqueN(rp_nj$WES_ID)
fwrite(rp_nj , "auto_candidate_analysis_proteomics10_11.csv")

rpsx <- rps[gene_class != "no rare variant" ]
rpsx <- rpsx[  ( SemSim > 40 | Semantic_sim > 2) ]
rpsx <- rpsx[ PROTEIN_ZSCORE < -2]

rpj <- rps[ WES_ID %like% "JAP" & ( SemSim > 40 | Semantic_sim > 2) ]
fwrite(rpj , "auto_candidate_analysis_proteomics10.csv")
exomiser <- fread("/data/agprokisch/working/dima/CowRD/data/processed_data/WES_WGS/exomiser/out/99302/exmsr_AD.variants.tsv")
