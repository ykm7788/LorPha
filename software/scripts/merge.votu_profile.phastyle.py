# -*- coding: utf-8 -*-
"""
Merge virus metadata and phase prediction TSV files by ContigID
Compatible with Python 3.11.4
Features:
- Robust error handling for file I/O and edge cases
- Memory-efficient data processing (no redundant storage)
- Clean column renaming and merging logic
- Preserves original order and handles missing values gracefully
"""

import argparse
import os
import csv
from typing import Dict, List, Optional


def read_virus_file(file_path: str) -> tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    Read virus metadata file (file1) and return headers + ContigID-indexed data
    Args:
        file_path: Path to virus TSV file
    Returns:
        Tuple of (header list, data dictionary {ContigID: row_data})
    """
    headers: List[str] = []
    data: Dict[str, Dict[str, str]] = {}
    
    # Handle file not found or permission errors
    try:
        with open(file_path, 'r', encoding='utf-8', newline='') as f:
            # Read TSV with DictReader for column name access
            reader = csv.DictReader(f, delimiter='\t')
            headers = reader.fieldnames or []
            
            # Skip if no headers or missing ContigID column
            if not headers or 'ContigID' not in headers:
                print(f"Warning: {file_path} has invalid headers (missing ContigID)", flush=True)
                return headers, data
            
            # Process rows (remove duplicates, strip whitespace)
            for row_num, row in enumerate(reader, start=2):  # Start counting at 2 (header=1)
                contig_id = row['ContigID'].strip()
                
                # Skip empty ContigID or duplicate entries
                if not contig_id or contig_id in data:
                    continue
                
                # Clean all values (strip whitespace)
                cleaned_row = {k: v.strip() for k, v in row.items()}
                data[contig_id] = cleaned_row
                
    except FileNotFoundError:
        print(f"Warning: Virus file {file_path} not found", flush=True)
    except PermissionError:
        print(f"Error: No permission to read {file_path}", flush=True)
    except Exception as e:
        print(f"Warning: Error reading {file_path} (row {row_num if 'row_num' in locals() else 'unknown'}): {str(e)}", flush=True)
    
    return headers, data


def process_phase_file(file_path: str) -> tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    Process phase prediction file (file2):
    1. Remove sequence_id column
    2. Rename columns to ContigID/VirusLife_PhaStyle/ScoreTemp_PhaStyle/ScoreViru_PhaStyle
    3. Return headers + ContigID-indexed data
    Args:
        file_path: Path to phase TSV file
    Returns:
        Tuple of (new header list, processed data dictionary {ContigID: row_data})
    """
    # Define new column names (after removing sequence_id)
    new_headers: List[str] = ['ContigID', 'VirusLife_PhaStyle', 'ScoreTemp_PhaStyle', 'ScoreViru_PhaStyle']
    data: Dict[str, Dict[str, str]] = {}
    
    # Original to new column mapping
    col_mapping = {
        'fasta_id': 'ContigID',
        'predicted_label': 'VirusLife_PhaStyle',
        'score_temperate': 'ScoreTemp_PhaStyle',
        'score_virulent': 'ScoreViru_PhaStyle'
    }
    
    try:
        with open(file_path, 'r', encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            orig_headers = reader.fieldnames or []
            
            # Validate required columns exist
            required_cols = list(col_mapping.keys())
            missing_cols = [col for col in required_cols if col not in orig_headers]
            if missing_cols:
                print(f"Warning: {file_path} missing required columns: {', '.join(missing_cols)}", flush=True)
                return new_headers, data
            
            # Process each row (skip duplicates, clean values)
            for row_num, row in enumerate(reader, start=2):
                # Extract and clean ContigID (fasta_id)
                contig_id = row['fasta_id'].strip()
                if not contig_id or contig_id in data:
                    continue
                
                # Map original columns to new names
                processed_row = {}
                for orig_col, new_col in col_mapping.items():
                    processed_row[new_col] = row[orig_col].strip()
                
                data[contig_id] = processed_row
                
    except FileNotFoundError:
        print(f"Warning: Phase file {file_path} not found", flush=True)
    except PermissionError:
        print(f"Error: No permission to read {file_path}", flush=True)
    except Exception as e:
        print(f"Warning: Error reading {file_path} (row {row_num if 'row_num' in locals() else 'unknown'}): {str(e)}", flush=True)
    
    return new_headers, data


def merge_datasets(
    virus_headers: List[str],
    virus_data: Dict[str, Dict[str, str]],
    phase_headers: List[str],
    phase_data: Dict[str, Dict[str, str]]
) -> tuple[List[str], List[Dict[str, str]]]:
    """
    Merge virus and phase datasets by ContigID (virus data as primary)
    - Preserves original order of virus data
    - Fills empty strings for missing phase data
    - Removes duplicate ContigID entries
    Args:
        virus_headers: Headers from virus file
        virus_data: ContigID-indexed virus data
        phase_headers: Processed phase headers
        phase_data: ContigID-indexed phase data
    Returns:
        Tuple of (merged headers, merged rows list)
    """
    # Combine headers (avoid duplicates)
    merged_headers = virus_headers + [h for h in phase_headers if h not in virus_headers]
    merged_rows: List[Dict[str, str]] = []
    
    # Merge rows (use virus data as base)
    for contig_id, virus_row in virus_data.items():
        merged_row = virus_row.copy()
        
        # Add phase data (empty string if missing)
        phase_row = phase_data.get(contig_id, {})
        for header in phase_headers:
            merged_row[header] = phase_row.get(header, '')
        
        merged_rows.append(merged_row)
    
    return merged_headers, merged_rows


def write_merged_file(output_dir: str, headers: List[str], rows: List[Dict[str, str]]) -> None:
    """
    Write merged data to TSV file (votu.meta.raw.txt)
    - Creates output directory if it doesn't exist
    - Handles permission errors gracefully
    Args:
        output_dir: Directory path for output file
        headers: Merged column headers
        rows: Merged data rows
    """
    # Create output directory if missing
    try:
        os.makedirs(output_dir, exist_ok=True)
    except PermissionError:
        print(f"Error: No permission to create directory {output_dir}", flush=True)
        return
    
    # Define full output path
    output_path = os.path.join(output_dir, 'votu.meta.raw.txt')
    
    # Write TSV file
    try:
        with open(output_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        
        print(f"Successfully wrote merged data to: {output_path}", flush=True)
        
        # Print summary stats
        print(f"Summary: {len(rows)} unique ContigID entries merged", flush=True)
        
    except PermissionError:
        print(f"Error: No permission to write to {output_path}", flush=True)
    except Exception as e:
        print(f"Error writing merged file: {str(e)}", flush=True)


def main() -> None:
    """Main workflow: parse args → read files → process → merge → write output"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Merge virus metadata and phase prediction TSV files by ContigID',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-v', '--virus', required=True, 
                        help='Path to virus metadata file (file1, TSV format with ContigID column)')
    parser.add_argument('-p', '--phase', required=True, 
                        help='Path to phase prediction file (file2, TSV format with fasta_id column)')
    parser.add_argument('-o', '--output', required=True, 
                        help='Output directory path (file will be saved as votu.meta.raw.txt)')
    
    args = parser.parse_args()

    # Step 1: Read and validate virus file
    virus_headers, virus_data = read_virus_file(args.virus)
    
    # Step 2: Process phase file (rename columns, remove sequence_id)
    phase_headers, phase_data = process_phase_file(args.phase)
    
    # Step 3: Merge datasets
    merged_headers, merged_rows = merge_datasets(virus_headers, virus_data, phase_headers, phase_data)
    
    # Step 4: Write merged output
    if merged_headers and merged_rows:
        write_merged_file(args.output, merged_headers, merged_rows)
    else:
        print("Warning: No valid data to merge - creating empty output file", flush=True)
        write_merged_file(args.output, merged_headers, [])


if __name__ == '__main__':
    main()