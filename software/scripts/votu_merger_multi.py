
import argparse
import os
import csv
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description='Merge votu coverage files')
    parser.add_argument('-d', required=True, help='Input directory')
    parser.add_argument('-o', required=True, help='Output directory')
    return parser.parse_args()


def process_file(filepath):
    data = defaultdict(dict)
    headers = []
    
    with open(filepath, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        
        for row in reader:
            if len(row) < 2:
                continue
            contig = row[0]
            for i in range(1, len(row)):
                try:
                    data[i][contig] = float(row[i]) if row[i] else 0
                except ValueError:
                    data[i][contig] = 0
    
    return headers, data


def main():
    args = parse_args()
    os.makedirs(args.o, exist_ok=True)
    
    files = sorted(Path(args.d).glob('*.votu.coverm.txt'))
    all_data = {}
    all_contigs = set()
    col_names = {}
    
    for file in files:
        sample = file.name.replace('.votu.coverm.txt', '')
        headers, data = process_file(file)
        all_data[sample] = data
        
        if not col_names:
            for i in data.keys():
                parts = headers[i].split(' ', 1)
                col_names[i] = parts[1].replace(' ', '_') if len(parts) > 1 else headers[i]
        
        for col_idx, contig_vals in data.items():
            all_contigs.update(contig_vals.keys())
    
    all_contigs = sorted(all_contigs)
    samples = sorted(all_data.keys())
    
    for col_idx, col_name in col_names.items():
        output_path = Path(args.o) / f'{col_name}.votu.txt'
        
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['Contig'] + samples)
            
            for contig in all_contigs:
                row = [contig]
                for sample in samples:
                    val = all_data[sample].get(col_idx, {}).get(contig, 0)
                    row.append(f'{val:.6f}')
                writer.writerow(row)


if __name__ == '__main__':
    main()