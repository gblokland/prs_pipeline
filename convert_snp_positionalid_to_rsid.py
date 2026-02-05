import requests
import csv
import argparse

def snp_to_rsid(chromosome, position, species='human'):
    server = "https://rest.ensembl.org"
    ext = f"/overlap/region/{species}/{chromosome}:{position}-{position}?feature=variation"
    headers = {"Content-Type": "application/json"}
    
    response = requests.get(server + ext, headers=headers)
    if not response.ok:
        return None

    data = response.json()
    if data:
        return [entry['id'] for entry in data if 'id' in entry]
    return None

def process_gwas_file(input_file, output_file, species='human'):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        header = next(reader)
        writer.writerow(['rsID'] + header[1:])  # Drop original SNP column

        for row in reader:
            if not row:
                continue

            snp_name = row[0]
            other_data = row[1:]

            try:
                chrom, pos, _ = snp_name.split(":")
                rsids = snp_to_rsid(chrom, pos, species=species)
                rsid_str = ", ".join(rsids) if rsids else "NA"
                writer.writerow([rsid_str] + other_data)
            except Exception as e:
                writer.writerow([f"Error: {e}"] + other_data)

def main():
    parser = argparse.ArgumentParser(description="Convert positional SNP IDs to rsIDs in a GWAS summary statistics file.")
    parser.add_argument("input_file", help="Path to the input GWAS summary statistics file (tab-delimited)")
    parser.add_argument("output_file", help="Path to the output file with rsIDs")
    parser.add_argument("--species", default="human", help="Species (default: human)")

    args = parser.parse_args()
    process_gwas_file(args.input_file, args.output_file, species=args.species)

if __name__ == "__main__":
    main()
