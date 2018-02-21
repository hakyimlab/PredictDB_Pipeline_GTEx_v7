import sys

def split_snp_annot(annot_file, out_prefix):
    snps_by_chr_files = [out_prefix + '.chr' + str(i) + '.txt' for i in range(1,23)]
    snps_by_chr = [open(f, 'w') for f in snps_by_chr_files]
    with open(annot_file, 'r') as ann:
        header = ann.readline().strip()
        for f in snps_by_chr:
            f.write(header + '\n')
        for line in ann:
            fields = line.strip().split()
            try:
                chrom = int(fields[0])
            except ValueError:
                continue
            snps_by_chr[chrom - 1].write(line.strip() + '\n')
    for f in snps_by_chr:
        f.close()

if __name__ == '__main__':
    annot_file = sys.argv[1]
    out_prefix = sys.argv[2]
    split_snp_annot(annot_file, out_prefix)  
