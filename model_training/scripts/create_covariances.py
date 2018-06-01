#!/bash/python

import gzip
import sqlite3
import pandas

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
             'Cells_Cultured_fibroblasts',
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

#tissues = ['Heart_Atrial_Appendage']

for tiss in tissues:
    try:
        print(tiss)
        # Get set of genes with significant models
        conn = sqlite3.connect('../dbs/gtex_v8_{}_qdir_signif.db'.format(tiss))
        c = conn.cursor()
        c.execute('select gene from extra')
        signif_genes = {row[0] for row in c.fetchall()}
        conn.close()
        with gzip.open('../covariances/gtex_v8_{}_qdir_signif.txt.gz'.format(tiss), 'w') as cov_out:
            hdr = ' '.join(['GENE', 'RSID1', 'RSID2', 'VALUE']) + '\n'
            cov_out.write(hdr.encode('utf-8'))
            for chrom in range(1,23):
                    c = pandas.read_table('../covariances/{}_nested_cv_chr{}_covariances.txt'.format(tiss, chrom), sep="\s+").drop_duplicates()
                    for row in c.itertuples():
                        gene = row.GENE
                        if gene in signif_genes:
                            l='{}\n'.format("\t".join([gene, row.RSID1, row.RSID2, str(row.VALUE)])).encode('utf-8')
                            cov_out.write(l)
    except:
        print("Error {}/{}".format(tiss, chrom))
