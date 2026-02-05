import pysam
import argparse
import csv

# Map chromosome number to RefSeq name
CHROM_MAP = {
    "1": "NC_000001.10", "2": "NC_000002.11", "3": "NC_000003.11",
    "4": "NC_000004.11", "5": "NC_000005.9",  "6": "NC_000006.11",
    "7": "NC_000007.13", "8": "NC_000008.10", "9": "NC_000009.11",
    "10": "NC_000010.10", "11": "NC_000011.9", "12": "NC_000012.11",
    "13": "NC_000013.10", "14": "NC_000014.8", "15": "NC_000015.9",
    "16": "NC_000016.9", "17": "NC_000017.10", "18": "NC_000018.9",
    "19": "NC_000019.9", "20": "NC_000020.10", "21": "NC_000021.8",
    "22": "NC_000022.10", "X": "NC_000023.10", "Y": "NC_000024.9"
}

def get_rsid(tabix_file, chrom, pos, ref, alt):
    mapped_chrom = CHROM_MAP.get(chrom)
    if not mapped_chrom:
        return "NA"

    try:
        records = tabix_file.fetch(mapped_chrom, int(pos) - 1, int(pos))
        for record in records:
            fields = record.strip().split('\t')
            vcf_pos = fields[1]
            rsid = fields[2]
            vcf_ref = fields[3]
            vcf_alt = fields[4].split(',')
            if vcf_pos == pos and vcf_ref == ref and alt in vcf_alt:
                return rsid
    except Exception as e:
        print(f"Error fetching {mapped_chrom}:{pos} - {e}")
    return "NA"

def convert_gwas(input_file, output_file, vcf_path):
    with pysam.TabixFile(vcf_path) as tabix, \
         open(input_file, 'r') as fin, \
         open(output_file, 'w', newline='') as fout:

        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')

        header = next(reader)
        writer.writerow(["rsid"] + header[1:])  # Replace SNP column with rsid

        for line in reader:
            snp = line[0]
            try:
                chrom, pos, alleles = snp.split(':')
                ref, alt = alleles.split('_')
            except ValueError:
                print(f"Skipping malformed SNP ID: {snp}")
                writer.writerow(["NA"] + line[1:])
                continue

            rsid = get_rsid(tabix, chrom, pos, ref, alt)
            writer.writerow([rsid] + line[1:])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert positional SNPs to rsIDs using a combined VCF.")
    parser.add_argument("--input", required=True, help="Input GWAS file (tab-separated, SNP in col 1)")
    parser.add_argument("--output", required=True, help="Output file with rsIDs")
    parser.add_argument("--vcf", required=True, help="Path to single bgzipped and indexed VCF file")
    args = parser.parse_args()

    convert_gwas(args.input, args.output, args.vcf)
