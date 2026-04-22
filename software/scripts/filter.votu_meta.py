# -*- coding: utf-8 -*-
"""
Virus Taxonomy Data Processor
Processes viral taxonomy data to generate standardized TSV output with priority-based filling.
Usage: python3 virus_tax_processor.py -i input.tsv -o output_directory
"""

import argparse
import os
import csv
import re
from typing import Dict, List


VIRUS_TAXONOMY_RANKS = [
    'superkingdom', 'clade', 'kingdom', 'phylum', 
    'class', 'order', 'family', 'genus', 'species'
]

INVALID_VALUES = {'', '-', 'Not found'}


def normalize_vhost_bacteria(tax_string: str) -> str:
    if 'Bacteria' not in tax_string:
        return tax_string

    tax_string = tax_string.replace('cellular organisms; Bacteria', 'd__Bacteria')

    parts = tax_string.split(';')
    if len(parts) < 2:
        return tax_string

    new_parts = [parts[0]]
    rank_prefixes = ['p__', 'c__', 'o__', 'f__', 'g__', 's__']

    for i, part in enumerate(parts[1:]):
        if i < len(rank_prefixes):
            clean_part = part.strip()
            if clean_part:
                new_parts.append(f"{rank_prefixes[i]}{clean_part}")
        else:
            new_parts.append(part.strip())

    result = ';'.join(new_parts)
    result = re.sub(r'__\s+', '__', result)

    return result


def parse_phabox_taxonomy(tax_string: str) -> Dict[str, str]:
    tax_data = {rank: '' for rank in VIRUS_TAXONOMY_RANKS}

    if not tax_string or tax_string.strip() in INVALID_VALUES:
        return tax_data

    for item in tax_string.strip().split(';'):
        item = item.strip()
        if not item or ':' not in item:
            continue

        rank_key, rank_value = item.split(':', 1)
        rank_key = rank_key.strip().lower()
        rank_value = rank_value.strip()

        if rank_key in tax_data and rank_value not in INVALID_VALUES:
            tax_data[rank_key] = rank_value

    return tax_data


def parse_list_taxonomy(tax_string: str) -> List[str]:
    if not tax_string or tax_string.strip() in INVALID_VALUES:
        return [''] * len(VIRUS_TAXONOMY_RANKS)

    tax_parts = [part.strip() for part in tax_string.strip().split(';')]

    if len(tax_parts) < len(VIRUS_TAXONOMY_RANKS):
        tax_parts += [''] * (len(VIRUS_TAXONOMY_RANKS) - len(tax_parts))

    return tax_parts[:len(VIRUS_TAXONOMY_RANKS)]


def fill_taxonomy_hierarchy(
    phabox_data: Dict[str, str],
    genomad_data: List[str],
    vhost_data: List[str]
) -> Dict[str, str]:
    final_taxonomy = phabox_data.copy()

    for idx, rank in enumerate(VIRUS_TAXONOMY_RANKS):
        if not final_taxonomy[rank] and genomad_data[idx] not in INVALID_VALUES:
            final_taxonomy[rank] = genomad_data[idx]

    for idx, rank in enumerate(VIRUS_TAXONOMY_RANKS):
        if not final_taxonomy[rank] and vhost_data[idx] not in INVALID_VALUES:
            final_taxonomy[rank] = vhost_data[idx]

    for rank in VIRUS_TAXONOMY_RANKS:
        if not final_taxonomy[rank]:
            final_taxonomy[rank] = f"unclassified_{rank}"

    return final_taxonomy


def extract_host_taxonomy(row: Dict[str, str]) -> str:
    gtdb_value = row.get('HostTax_phaboxGTDB', '').strip()
    if gtdb_value and gtdb_value not in INVALID_VALUES:
        return gtdb_value

    ncbi_value = row.get('HostTax_phaboxNCBI', '').strip()
    if ncbi_value and ncbi_value not in {'', '-'}:
        return ncbi_value

    vhost_value = row.get('HostTax_vhost', '').strip()
    return vhost_value if vhost_value else ''


def determine_virus_lifestyle(row: Dict[str, str]) -> str:
    try:
        score_temp = float(row.get('ScoreTemp_PhaStyle', 0.0)) if row.get('ScoreTemp_PhaStyle', '').strip() else 0.0
        score_viru = float(row.get('ScoreViru_PhaStyle', 0.0)) if row.get('ScoreViru_PhaStyle', '').strip() else 0.0

        if score_temp > 0.7:
            return "temperate"
        elif score_viru > 0.7:
            return "virulent"

    except (ValueError, TypeError):
        pass

    lifestyle_phabox = row.get('VirusLife_phabox', '').strip()
    return lifestyle_phabox if lifestyle_phabox not in INVALID_VALUES else ''


def validate_input_columns(reader: csv.DictReader) -> None:
    required_columns = [
        'ContigID', 'VirusTax_phabox', 'VirusTax_genomad', 'VirusTax_vhost',
        'HostTax_phaboxGTDB', 'HostTax_phaboxNCBI', 'HostTax_vhost',
        'VirusLife_phabox', 'ScoreTemp_PhaStyle', 'ScoreViru_PhaStyle'
    ]

    missing_columns = [col for col in required_columns if col not in reader.fieldnames]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")


def process_tsv(input_path: str, output_path: str) -> None:
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    processed_rows = []

    try:
        with open(input_path, 'r', encoding='utf-8', newline='') as infile:
            reader = csv.DictReader(infile, delimiter='\t')
            validate_input_columns(reader)

            for row_num, row in enumerate(reader, start=2):
                contig_id = row['ContigID'].strip()

                if not contig_id:
                    continue

                if 'HostTax_vhost' in row and 'Bacteria' in row['HostTax_vhost']:
                    row['HostTax_vhost'] = normalize_vhost_bacteria(row['HostTax_vhost'])

                phabox_tax = parse_phabox_taxonomy(row['VirusTax_phabox'])
                genomad_tax = parse_list_taxonomy(row['VirusTax_genomad'])
                vhost_tax = parse_list_taxonomy(row['VirusTax_vhost'])

                final_tax = fill_taxonomy_hierarchy(phabox_tax, genomad_tax, vhost_tax)

                virus_taxonomy = ';'.join([final_tax[rank] for rank in VIRUS_TAXONOMY_RANKS])

                host_taxonomy = extract_host_taxonomy(row)
                virus_lifestyle = determine_virus_lifestyle(row)

                processed_rows.append({
                    'ContigID': contig_id,
                    'VirusTaxonomy': virus_taxonomy,
                    'HostTaxonomy': host_taxonomy,
                    'VirusLifestyle': virus_lifestyle
                })

        with open(output_path, 'w', encoding='utf-8', newline='') as outfile:
            output_fields = ['ContigID', 'VirusTaxonomy', 'HostTaxonomy', 'VirusLifestyle']
            writer = csv.DictWriter(outfile, fieldnames=output_fields, delimiter='\t')

            writer.writeheader()
            writer.writerows(processed_rows)

    except FileNotFoundError:
        raise FileNotFoundError(f"Input file not found: {input_path}")
    except Exception as e:
        raise RuntimeError(f"Error processing row {row_num if 'row_num' in locals() else 'unknown'}: {str(e)}")


def main():
    parser = argparse.ArgumentParser(
        description='Process viral taxonomy data to generate standardized outputs'
    )

    parser.add_argument('-i', '--input', 
                        required=True,
                        help='Path to input TSV file')

    parser.add_argument('-o', '--output', 
                        required=True,
                        help='Output directory path')

    args = parser.parse_args()
    output_file = os.path.join(args.output, 'votu.meta.filter.txt')

    try:
        process_tsv(args.input, output_file)
    except Exception as e:
        print(f"Error: {str(e)}", flush=True)
        exit(1)


if __name__ == '__main__':
    main()