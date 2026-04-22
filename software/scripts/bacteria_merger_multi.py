import argparse
import os
import csv
import re
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description='Process tax coverage files')
    parser.add_argument('-d', required=True, help='Input directory')
    parser.add_argument('-o', required=True, help='Output directory')
    return parser.parse_args()


def normalize_taxa(taxa):
    taxa = re.sub(r'^d_', 'd__', taxa)
    taxa = re.sub(r';p_', ';p__', taxa)
    taxa = re.sub(r';c_', ';c__', taxa)
    taxa = re.sub(r';o_', ';o__', taxa)
    taxa = re.sub(r';f_', ';f__', taxa)
    taxa = re.sub(r';g_', ';g__', taxa)
    taxa = re.sub(r';s_', ';s__', taxa)
    return taxa


def process_file(filepath):
    data = defaultdict(lambda: defaultdict(float))
    headers = []
    
    with open(filepath, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        
        for row in reader:
            if len(row) < 4:
                continue
            identity = float(row[1]) if row[1] else 0
            taxa = normalize_taxa(row[2]) if identity >= 0.8 else 'unclassified'
            
            for i in range(3, len(row)):
                try:
                    val = float(row[i]) if row[i] else 0
                    data[i][taxa] += val
                except ValueError:
                    continue
    
    return headers, data


def main():
    args = parse_args()
    os.makedirs(args.o, exist_ok=True)
    
    files = sorted(Path(args.d).glob('*.tax.coverm.txt'))
    all_data = {}
    all_taxa = set()
    col_headers = {}
    
    for file in files:
        sample = file.name.replace('.tax.coverm.txt', '')
        headers, data = process_file(file)
        all_data[sample] = data
        
        if not col_headers:
            col_headers = {i: headers[i].replace(' ', '_') for i in data.keys()}
        
        for col_idx, taxa_vals in data.items():
            all_taxa.update(taxa_vals.keys())
    
    all_taxa = sorted(all_taxa)
    samples = sorted(all_data.keys())
    
    for col_idx, col_name in col_headers.items():
        output_path = Path(args.o) / f'{col_name}.prokaryote.txt'
        
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['taxa'] + samples)
            
            for taxa in all_taxa:
                row = [taxa]
                for sample in samples:
                    val = all_data[sample].get(col_idx, {}).get(taxa, 0)
                    row.append(f'{val:.6f}' if isinstance(val, float) else str(val))
                writer.writerow(row)


if __name__ == '__main__':
    main()