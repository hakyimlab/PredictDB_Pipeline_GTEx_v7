library(dplyr)
library(RSQLite)
"%&%" <- function(a,b) paste(a,b, sep='')

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
driver <- dbDriver('SQLite')
gene_annot <- read.table("../../prepare_data/expression/gencode.v19.genes.patched_contigs.parsed.txt", header = T, stringsAsFactors = F)

for (tiss in tissues) {
  print(tiss)
  # Extra table ----
  model_summaries <- read.table('../summary/' %&% tiss %&% '_nested_cv_chr1_model_summaries.txt', header = T, stringsAsFactors = F)
  tiss_summary <- read.table('../summary/' %&% tiss %&% '_nested_cv_chr1_tiss_chr_summary.txt', header = T, stringsAsFactors = F)
  
  n_samples <- tiss_summary$n_samples
  
  for (i in 2:22) {
    model_summaries <- rbind(model_summaries,
                             read.table('../summary/' %&% tiss %&% '_nested_cv_chr' %&% as.character(i) %&% '_model_summaries.txt', header = T, stringsAsFactors = F))
    tiss_summary <- rbind(tiss_summary,
                             read.table('../summary/' %&% tiss %&% '_nested_cv_chr' %&% as.character(i) %&% '_tiss_chr_summary.txt', header = T, stringsAsFactors = F))
  }
  
  model_summaries <- rename(model_summaries, gene = gene_id)

  conn <- dbConnect(drv = driver, '../dbs/gtex_v7_' %&% tiss %&% '_imputed_europeans_tw0.5.db')
  dbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)
  dbGetQuery(conn, "CREATE INDEX gene_model_summary ON model_summaries (gene)")
  
  # Weights Table -----
  weights <- read.table('../weights/' %&% tiss %&% '_nested_cv_chr1_weights.txt', header = T, stringsAsFactors = F)
  for (i in 2:22) {
    weights <- rbind(weights,
                       read.table('../weights/' %&% tiss %&% '_nested_cv_chr' %&% as.character(i) %&% '_weights.txt', header = T, stringsAsFactors = F))
  }
  weights <- rename(weights, gene = gene_id)
  dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
  dbGetQuery(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbGetQuery(conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbGetQuery(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  
  # Sample_info Table ----
  sample_info <- data.frame(n_samples = n_samples, population = 'europeans', tissue = tiss)
  dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)
  
  # Construction Table ----
  construction <- tiss_summary %>%
                    select(chrom, cv_seed) %>%
                    rename(chromosome = chrom)
  dbWriteTable(conn, 'construction', construction, overwrite = TRUE)
}
