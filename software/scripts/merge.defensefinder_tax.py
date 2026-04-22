import argparse
import os
import csv
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description='Merge defense finder and taxonomy results')
    parser.add_argument('-d', required=True, help='Directory containing defense finder .tsv files')
    parser.add_argument('-t', required=True, help='Directory containing taxonomy .txt files')
    parser.add_argument('-o', required=True, help='Output directory')
    return parser.parse_args()


def get_sample_files(directory, suffix):
    files = {}
    for f in Path(directory).glob(f'*{suffix}'):
        sample = f.name.replace(suffix, '')
        files[sample] = f
    return files


def normalize_tax(tax):
    tax = re.sub(r'^d_', 'd__', tax)
    tax = re.sub(r';p_', ';p__', tax)
    tax = re.sub(r';c_', ';c__', tax)
    tax = re.sub(r';o_', ';o__', tax)
    tax = re.sub(r';f_', ';f__', tax)
    tax = re.sub(r';g_', ';g__', tax)
    tax = re.sub(r';s_', ';s__', tax)
    return tax


def load_tax_file(filepath):
    data = {}
    with open(filepath, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) >= 9:
                data[row[0]] = (row[7], normalize_tax(row[8]))
    return data


def main():
    args = parse_args()
    os.makedirs(args.o, exist_ok=True)

    df_files = get_sample_files(args.d, '.cluster_defense_finder_genes.tsv')
    tax_files = get_sample_files(args.t, '.tax.txt')
    common_samples = sorted(set(df_files) & set(tax_files))

    output_path = Path(args.o) / 'Defense_finder.tax.txt'
    with open(output_path, 'w', newline='', encoding='utf-8') as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerow(['sample', 'ContigID', 'GeneID', 'DF_gene', 'type', 'subtype', 'Tax_identity', 'Tax'])

        for sample in common_samples:
            tax_data = load_tax_file(tax_files[sample])
            with open(df_files[sample], 'r', newline='', encoding='utf-8') as inf:
                reader = csv.DictReader(inf, delimiter='\t')
                for row in reader:
                    replicon = row['replicon']
                    tax_id, tax = tax_data.get(replicon, ('', ''))
                    writer.writerow([sample, replicon, row['hit_id'], row['gene_name'],
                                   row['type'], row['subtype'], tax_id, tax])


if __name__ == '__main__':
    main()