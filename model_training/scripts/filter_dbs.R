library(dplyr)
library(RSQLite)
"%&%" <- function(a,b) paste(a,b,sep='')

driver <- dbDriver("SQLite")

tissues <- c('Adipose_Subcutaneous',
             'Adipose_Visceral_Omentum',
             'Adrenal_Gland',
             'Artery_Aorta',
             'Artery_Coronary',
             'Artery_Tibial',
             'Brain_Amygdala',
             'Brain_Anterior_cingulate_cortex_BA24',
             'Brain_Caudate_basal_ganglia',
             'Brain_Cerebellar_Hemisphere',
             'Brain_Cerebellum',
             'Brain_Cortex',
             'Brain_Frontal_Cortex_BA9',
             'Brain_Hippocampus',
             'Brain_Hypothalamus',
             'Brain_Nucleus_accumbens_basal_ganglia',
             'Brain_Putamen_basal_ganglia',
             'Brain_Spinal_cord_cervical_c-1',
             'Brain_Substantia_nigra',
             'Breast_Mammary_Tissue',
             'Cells_EBV-transformed_lymphocytes',
             'Cells_Transformed_fibroblasts',
             'Colon_Sigmoid',
             'Colon_Transverse',
             'Esophagus_Gastroesophageal_Junction',
             'Esophagus_Mucosa',
             'Esophagus_Muscularis',
             'Heart_Atrial_Appendage',
             'Heart_Left_Ventricle',
             'Liver',
             'Lung',
             'Minor_Salivary_Gland',
             'Muscle_Skeletal',
             'Nerve_Tibial',
             'Ovary',
             'Pancreas',
             'Pituitary',
             'Prostate',
             'Skin_Not_Sun_Exposed_Suprapubic',
             'Skin_Sun_Exposed_Lower_leg',
             'Small_Intestine_Terminal_Ileum',
             'Spleen',
             'Stomach',
             'Testis',
             'Thyroid',
             'Uterus',
             'Vagina',
             'Whole_Blood'
)

for (tissue in tissues) {
  print(tissue)
  unfiltered_db <- "../dbs/gtex_v7_" %&% tissue %&% "_imputed_europeans_tw0.5.db"
  filtered_db <- "../dbs/gtex_v7_" %&% tissue %&% "_imputed_europeans_tw_0.5_signif.db"
  in_conn <- dbConnect(driver, unfiltered_db)
  out_conn <- dbConnect(driver, filtered_db)
  model_summaries <- dbGetQuery(in_conn, 'select * from model_summaries where zscore_pval < 0.05 and rho_avg > 0.1')
  model_summaries <- model_summaries %>% rename(pred.perf.R2 = rho_avg_squared, genename = gene_name, pred.perf.pval = zscore_pval, n.snps.in.model = n_snps_in_model)
  model_summaries$pred.perf.qval <- NA
  dbWriteTable(out_conn, 'extra', model_summaries)
  construction <- dbGetQuery(in_conn, 'select * from construction')
  dbWriteTable(out_conn, 'construction', construction)
  sample_info <- dbGetQuery(in_conn, 'select * from sample_info')
  dbWriteTable(out_conn, 'sample_info', sample_info)
  weights <- dbGetQuery(in_conn, 'select * from weights')
  weights <- weights %>% filter(gene %in% model_summaries$gene) %>% rename(eff_allele = alt, ref_allele = ref, weight = beta)
  dbWriteTable(out_conn, 'weights', weights)
  dbGetQuery(out_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbGetQuery(out_conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbGetQuery(out_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbGetQuery(out_conn, "CREATE INDEX gene_model_summary ON extra (gene)")
}
