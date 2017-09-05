from collections import defaultdict
import gzip

def get_ids_from_vcf(vcf_file):
    with gzip.open(vcf_file) as fh:
        for line in fh:
            line = line.decode('utf-8')
            if not line.startswith("##"):
                break
        sample_ids = line.strip().split()[9:]
        return set(sample_ids)

vcf_file = '../v7/56055/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2016-01-15_v7/genotypes/WGS/variant_calls/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz'
vcf_sample_ids = get_ids_from_vcf(vcf_file)
tissue_donors = defaultdict(list)

samples_file = '../v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SampleAttributesDS.txt'

with open(samples_file) as fh:
    hdr = fh.readline() # Drop header
    for line in fh:
        if line == '' or line == '\n':
            continue
        fields = line.strip().split('\t')
        sample_id = fields[0]
        tissue = fields[13]
        frz = fields[27]
        donor_id = '-'.join(sample_id.split('-')[:2])
        if donor_id in vcf_sample_ids and frz == 'RNASEQ':
            tissue_donors[tissue].append(sample_id + '\n')

for tissue in tissue_donors:
    with open(tissue.replace(' - ', '_').replace(' ', '_') + '_donors.txt', 'w') as fh:
        fh.writelines(tissue_donors[tissue])

