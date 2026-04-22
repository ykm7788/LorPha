import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(description='Merge taxonomy and coverage data')
    parser.add_argument('-t', required=True, help='Taxonomy file (tab-delimited, no header)')
    parser.add_argument('-c', required=True, help='Coverage file (tab-delimited, with header)')
    parser.add_argument('-o', required=True, help='Output file')
    return parser.parse_args()


def load_tax(filepath):
    data = {}
    with open(filepath, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) >= 10:
                data[row[0]] = (row[8], row[9])
    return data


def process_header(header):
    new_header = [header[0]]
    for col in header[1:]:
        parts = col.split(' ', 1)
        new_header.append(parts[1] if len(parts) > 1 else col)
    return new_header


def main():
    args = parse_args()
    
    tax_data = load_tax(args.t)
    
    with open(args.c, 'r', newline='', encoding='utf-8') as cf, \
         open(args.o, 'w', newline='', encoding='utf-8') as outf:
        reader = csv.reader(cf, delimiter='\t')
        writer = csv.writer(outf, delimiter='\t')
        
        header = next(reader)
        new_header = process_header(header)
        writer.writerow(['Contig', 'identity', 'taxa'] + new_header[1:])
        
        for row in reader:
            contig = row[0]
            identity, taxa = tax_data.get(contig, ('', ''))
            writer.writerow([contig, identity, taxa] + row[1:])


if __name__ == '__main__':
    main()

