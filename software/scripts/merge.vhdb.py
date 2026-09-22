import argparse
import csv
import os

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--votu', required=True, help='path to votu.meta.filter file')
    parser.add_argument('-d', '--db', required=True, help='path to virushostdb.tsv file')
    args = parser.parse_args()

    # Build a lookup from virus name to (host lineage, host name)
    virus2host = {}
    with open(args.db, encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        v_idx = header.index('virus name')
        n_idx = header.index('host name')
        l_idx = header.index('host lineage')
        for row in reader:
            if len(row) <= max(v_idx, n_idx, l_idx):
                continue
            key = row[v_idx].strip()
            virus2host[key] = (row[l_idx].strip(), row[n_idx].strip())

    # Read the votu file
    with open(args.votu, encoding='utf-8') as f:
        lines = f.read().splitlines()

    header = lines[0].split('\t')
    ht_idx = header.index('HostTaxonomy')
    vt_idx = header.index('VirusTaxonomy')

    out_lines = [lines[0]]
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split('\t')
        if parts[ht_idx].strip() == '':
            # Extract the text after the last ';' in VirusTaxonomy
            vt = parts[vt_idx].strip()
            key = vt.split(';')[-1].strip()
            if key in virus2host:
                lineage, host_name = virus2host[key]
                # Split host lineage on ';' and strip whitespace
                elems = [e.strip() for e in lineage.split(';')]

                # Map 1-based positions 2,4,5,6,7,8 to domain/phylum/class/order/family/genus
                def get(i):
                    return elems[i - 1] if len(elems) >= i else 'unclassified'

                taxa = ';'.join([
                    'd__' + get(2),
                    'p__' + get(4),
                    'c__' + get(5),
                    'o__' + get(6),
                    'f__' + get(7),
                    'g__' + get(8),
                    's__' + host_name,
                ])
                parts[ht_idx] = taxa
                line = '\t'.join(parts)
        out_lines.append(line)

    out_path = os.path.join(os.path.dirname(args.votu), 'votu.meta.taxa.txt')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.  join(out_lines) + '\n')

    print(out_path)

if __name__ == '__main__':
    main()
