#!/bash/python

import gzip
import sqlite3

tissues = ['Adipose_Subcutaneous',
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
             'Whole_Blood']

for tiss in tissues:
    print(tiss)
    # Get set of genes with significant models
    conn = sqlite3.connect('../dbs/gtex_v7_{}_imputed_europeans_tw_0.5_signif.db'.format(tiss))
    c = conn.cursor()
    c.execute('select gene from extra')
    signif_genes = {row[0] for row in c.fetchall()}
    conn.close()
    with gzip.open('../covariances/gtex_v7_{}_imputed_eur_covariances.txt.gz'.format(tiss), 'w') as cov_out:
        hdr = ' '.join(['GENE', 'RSID1', 'RSID2', 'VALUE']) + '\n'
        cov_out.write(hdr.encode('utf-8'))
        for chrom in range(1,23):
            with open('../covariances/{}_nested_cv_chr{}_covariances.txt'.format(tiss, chrom)) as cov_in:
                # Skip header
                line = cov_in.readline()
                for line in cov_in:
                    gene = line.strip().split()[0]
                    if gene in signif_genes:
                        cov_out.write((line.strip() +'\n').encode('utf-8')) 

