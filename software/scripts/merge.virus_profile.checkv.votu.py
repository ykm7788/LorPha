import pandas as pd
import argparse
import os
import sys

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Merge three files based on contig IDs from file3')
    parser.add_argument('-p', '--file1', required=True, help='Path to virus.txt (tab-separated with header)')
    parser.add_argument('-c', '--file2', required=True, help='Path to checkv (tab-separated with header)')
    parser.add_argument('-v', '--file3', required=True, help='Path to votu.contig (no header, one contig ID per line)')
    parser.add_argument('-o', '--output', required=True, help='Output directory path')
    args = parser.parse_args()

    # -------------------------- 1. Read and validate file1 --------------------------
    def read_file1(file_path):
        """Read file1 with error handling and validation"""
        try:
            # Read tab-separated file, keep empty strings instead of NaN
            df = pd.read_csv(file_path, sep='\t', dtype=str, keep_default_na=False)
            # Check if file is empty or only has header
            if df.empty or len(df.index) == 0:
                print(f"Warning: {file_path} is empty or only contains header", file=sys.stderr)
                sys.exit(1)
            return df
        except FileNotFoundError:
            print(f"Warning: {file_path} not found", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Warning: Error reading {file_path}: {str(e)}", file=sys.stderr)
            sys.exit(1)

    # -------------------------- 2. Read and validate file2 --------------------------
    def read_file2(file_path):
        """Read file2, keep specified columns with error handling"""
        keep_cols = ['contig_id', 'provirus', 'checkv_quality', 'miuvig_quality', 'contamination']
        try:
            df = pd.read_csv(file_path, sep='\t', dtype=str, keep_default_na=False)
            if df.empty or len(df.index) == 0:
                print(f"Warning: {file_path} is empty or only contains header", file=sys.stderr)
                sys.exit(1)
            # Add missing columns with empty values
            for col in keep_cols:
                if col not in df.columns:
                    df[col] = ''
            # Rename columns as required
            df = df[keep_cols].rename(columns={
                'provirus': 'checkv_provirus',
                'contamination': 'checkv_contamination'
            })
            return df
        except FileNotFoundError:
            print(f"Warning: {file_path} not found", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Warning: Error reading {file_path}: {str(e)}", file=sys.stderr)
            sys.exit(1)

    # -------------------------- 3. Read and validate file3 --------------------------
    def read_file3(file_path):
        """Read file3 (contig IDs) with error handling and deduplication"""
        try:
            # Read without header, skip blank lines
            df = pd.read_csv(
                file_path, sep='\t', header=None, names=['ContigID'], 
                dtype=str, skip_blank_lines=True, keep_default_na=False
            )
            # Remove empty ContigID rows and duplicates
            df = df[df['ContigID'] != ''].drop_duplicates(subset=['ContigID'], keep='first')
            if df.empty:
                print(f"Warning: {file_path} has no valid contig IDs", file=sys.stderr)
                sys.exit(1)
            return df
        except FileNotFoundError:
            print(f"Warning: {file_path} not found", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Warning: Error reading {file_path}: {str(e)}", file=sys.stderr)
            sys.exit(1)

    # Execute file reading with validation
    df1 = read_file1(args.file1)
    df2 = read_file2(args.file2)
    df3 = read_file3(args.file3)

    # -------------------------- 4. Merge all dataframes --------------------------
    # Merge file3 (base) with file1 (left join to keep only file3 IDs)
    merged_df = pd.merge(df3, df1, on='ContigID', how='left').fillna('')
    # Rename file2 key column and merge
    df2_renamed = df2.rename(columns={'contig_id': 'ContigID'})
    merged_df = pd.merge(merged_df, df2_renamed, on='ContigID', how='left').fillna('')
    # Final deduplication to remove duplicate rows
    merged_df = merged_df.drop_duplicates(subset=['ContigID'], keep='first')

    # -------------------------- 5. Output results --------------------------
    # Create output directory if not exists
    os.makedirs(args.output, exist_ok=True)
    output_file = os.path.join(args.output, 'votu.profile.txt')
    # Save tab-separated file without index, empty values as blank
    merged_df.to_csv(output_file, sep='\t', index=False, na_rep='')

if __name__ == "__main__":
    main()