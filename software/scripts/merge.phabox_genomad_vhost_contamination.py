"""
Processes three viral analysis output files (PhaBox, geNomad, vhost),
applies specific filters, merges them based on ContigID, and writes the
results to two output files.
"""

import pandas as pd
import argparse
import os

def check_file_status(file_path, file_name):
    """
    Checks if a file exists, is empty, or contains only headers.

    Args:
        file_path (str): Path to the file to check.
        file_name (str): A descriptive name for the file ("file1", "file2", etc.).

    Returns:
        str: The status of the file ("valid", "missing", "empty").
    """
    if not os.path.exists(file_path):
        print(f"Warning: {file_name} is missing")
        return "missing"
    if os.path.getsize(file_path) == 0:
        print(f"Warning: {file_name} is empty")
        return "empty"
    try:
        # Read a single row to check for content beyond headers
        sample_kwargs = {"sep": '\t', "dtype": str, "nrows": 1}
        if file_name == "file3":
             sample_kwargs["header"] = None
        sample = pd.read_csv(file_path, **sample_kwargs)
        if sample.empty:
            print(f"Warning: {file_name} has only headers or no data")
            return "empty"
    except pd.errors.EmptyDataError:
        print(f"Warning: {file_name} is empty")
        return "empty"
    return "valid"


def process_file1(file_path):
    """
    Processes the PhaBox output file (file1).

    Filters for rows where Pred (column 3) is 'virus' AND PhaMerConfidence (column 6) is 'high-confidence'.

    Args:
        file_path (str): Path to file1.

    Returns:
        pandas.DataFrame: Processed and filtered DataFrame with specific columns,
                          or an empty DataFrame with the correct columns if processing fails.
    """
    status = check_file_status(file_path, "file1")
    if status in ["missing", "empty"]:
        return pd.DataFrame(columns=[
            'ContigID', 'PhaMerConfidence', 'VirusTax_phabox',
            'VirusLife_phabox', 'HostTax_phaboxNCBI', 'HostTax_phaboxGTDB'
        ])

    try:
        df = pd.read_csv(file_path, sep='\t', dtype=str)
    except Exception as e:
        # Error during read is treated as failure, return empty DF
        return pd.DataFrame(columns=[
            'ContigID', 'PhaMerConfidence', 'VirusTax_phabox',
            'VirusLife_phabox', 'HostTax_phaboxNCBI', 'HostTax_phaboxGTDB'
        ])

    if df.empty or df.shape[1] < 17:
        print(f"Warning: file1 does not have enough columns or is empty after reading.")
        return pd.DataFrame(columns=[
            'ContigID', 'PhaMerConfidence', 'VirusTax_phabox',
            'VirusLife_phabox', 'HostTax_phaboxNCBI', 'HostTax_phaboxGTDB'
        ])

    # Identify rows where Pred column (index 2) is 'virus'
    pred_is_virus_mask = df.iloc[:, 2] == 'virus'

    # Get the relevant columns first (for easier subsequent filtering)
    # These correspond to the original 1st, 6th, 7th, 11th, 16th, 17th columns
    temp_df = df.loc[pred_is_virus_mask, df.columns[[0, 5, 6, 11, 16, 17]]].copy()
    temp_df.columns = [
        'ContigID', 'PhaMerConfidence', 'VirusTax_phabox',
        'VirusLife_phabox', 'HostTax_phaboxNCBI', 'HostTax_phaboxGTDB'
    ]

    if temp_df.empty:
        print("Warning: file1 yielded no valid data after initial 'virus' prediction filter.")
        return pd.DataFrame(columns=[
            'ContigID', 'PhaMerConfidence', 'VirusTax_phabox',
            'VirusLife_phabox', 'HostTax_phaboxNCBI', 'HostTax_phaboxGTDB'
        ])

    # Apply the filter: PhaMerConfidence must be 'high-confidence'
    final_filter_mask = temp_df['PhaMerConfidence'] == 'high-confidence'

    # Apply the final filter
    filtered_result_df = temp_df[final_filter_mask]

    if filtered_result_df.empty:
        print("Warning: file1 yielded no valid data after applying the high-confidence filtering criterion.")
    
    return filtered_result_df


def process_file2(file_path):
    """
    Processes the geNomad output file (file2).

    Simply selects the ContigID (1st) and VirusTax_genomad (11th) columns.

    Args:
        file_path (str): Path to file2.

    Returns:
        pandas.DataFrame: Processed DataFrame with selected columns,
                          or an empty DataFrame with the correct columns if processing fails.
    """
    status = check_file_status(file_path, "file2")
    if status in ["missing", "empty"]:
        return pd.DataFrame(columns=['ContigID', 'VirusTax_genomad'])

    try:
        df = pd.read_csv(file_path, sep='\t', dtype=str)
    except Exception as e:
        # Error during read is treated as failure, return empty DF
        return pd.DataFrame(columns=['ContigID', 'VirusTax_genomad'])

    if not df.empty:
        # Select ContigID (0) and VirusTax_genomad (10) columns
        selected_cols = df.iloc[:, [0, 10]].copy()
        selected_cols.columns = ['ContigID', 'VirusTax_genomad']

        if selected_cols.empty:
            print("Warning: file2 yielded no valid data after column selection.")
            return pd.DataFrame(columns=['ContigID', 'VirusTax_genomad'])
        
        return selected_cols
    else:
        print("Warning: file2 yielded no valid data after reading.")
        return pd.DataFrame(columns=['ContigID', 'VirusTax_genomad'])


def process_file3(file_path):
    """
    Processes the vhost output file (file3).

    Splits the 21st column by '|' and takes the first four parts for virus/host info.
    Selects the ContigID (1st) and the processed 21st column.

    Args:
        file_path (str): Path to file3.

    Returns:
        pandas.DataFrame: Processed DataFrame with ContigID and parsed virus/host/taxonomy columns,
                          or an empty DataFrame with the correct columns if processing fails.
    """
    status = check_file_status(file_path, "file3")
    if status in ["missing", "empty"]:
        return pd.DataFrame(columns=[
            'ContigID', 'Virus_vhost', 'Host_vhost', 'VirusTax_vhost', 'HostTax_vhost'
        ])

    try:
        df = pd.read_csv(file_path, sep='\t', dtype=str, header=None)
    except Exception as e:
        # Error during read is treated as failure, return empty DF
        return pd.DataFrame(columns=[
            'ContigID', 'Virus_vhost', 'Host_vhost', 'VirusTax_vhost', 'HostTax_vhost'
        ])

    if not df.empty and df.shape[1] > 20:
        # Select ContigID (0) and the raw column to be split (20)
        selected_data = df.iloc[:, [0, 20]].copy()

        # Prepare the column to be split: remove trailing '||' and split by '|'
        split_series = selected_data.iloc[:, 1].str.rstrip('||').str.split('|', expand=True)

        # Ensure the resulting DataFrame has exactly 4 columns
        for i in range(split_series.shape[1], 4):
            split_series[i] = '' # Pad with empty string if fewer than 4 parts
        split_series = split_series.iloc[:, :4] # Truncate if more than 4 parts
        split_series.columns = ['Virus_vhost', 'Host_vhost', 'VirusTax_vhost', 'HostTax_vhost']

        # Combine the ContigID column with the split results
        result = pd.concat([selected_data.iloc[:, [0]], split_series], axis=1)
        result.columns = ['ContigID', 'Virus_vhost', 'Host_vhost', 'VirusTax_vhost', 'HostTax_vhost']

        if result.empty:
            print("Warning: file3 yielded no valid data after parsing.")
            return pd.DataFrame(columns=[
                'ContigID', 'Virus_vhost', 'Host_vhost', 'VirusTax_vhost', 'HostTax_vhost'
            ])
        
        return result
    else:
        print("Warning: file3 yielded no valid data after reading or had insufficient columns.")
        return pd.DataFrame(columns=[
            'ContigID', 'Virus_vhost', 'Host_vhost', 'VirusTax_vhost', 'HostTax_vhost'
        ])


def main():
    """Main function to parse arguments, process files, merge results, and write outputs."""
    parser = argparse.ArgumentParser(description='Process viral data files (PhaBox, geNomad, vhost)')
    parser.add_argument('-p', '--file1', required=True, help='Path to file1 (PhaBox output)')
    parser.add_argument('-g', '--file2', required=True, help='Path to file2 (geNomad output)')
    parser.add_argument('-v', '--file3', required=True, help='Path to file3 (vhost output)')
    parser.add_argument('-o', '--output', required=True, help='Output directory path')
    args = parser.parse_args()

    # Process all input files
    df1 = process_file1(args.file1)
    df2 = process_file2(args.file2)
    df3 = process_file3(args.file3)

    # Start with df1's results
    merged_df = df1
    
    # Outer join with df2 on ContigID if df2 is not empty
    if not df2.empty:
        merged_df = pd.merge(merged_df, df2, on='ContigID', how='outer')

    # Outer join with df3 on ContigID if df3 is not empty
    if not df3.empty:
        merged_df = pd.merge(merged_df, df3, on='ContigID', how='outer')

    # Replace any remaining NaN values with empty strings
    merged_df = merged_df.fillna('')
    
    # Remove any duplicate rows that might have been introduced during merging
    merged_df = merged_df.drop_duplicates()

    # Define the required columns for the final output
    required_columns = [
        'ContigID', 'PhaMerConfidence', 'VirusTax_phabox', 'VirusLife_phabox',
        'HostTax_phaboxNCBI', 'HostTax_phaboxGTDB', 'VirusTax_genomad',
        'Virus_vhost', 'Host_vhost', 'VirusTax_vhost', 'HostTax_vhost'
    ]

    # Ensure all required columns exist in the final DataFrame.
    # If a column is missing, add it with empty string values.
    for col in required_columns:
        if col not in merged_df.columns:
            merged_df[col] = ''

    # Reorder the DataFrame columns to match the required order
    merged_df = merged_df[required_columns]

    # Create the output directory if it doesn't already exist
    os.makedirs(args.output, exist_ok=True)

    # Write the main merged output file
    main_output_path = os.path.join(args.output, 'virus.txt')
    merged_df.to_csv(main_output_path, sep='\t', index=False)

    # Write the ContigID-only output file
    contig_output_path = os.path.join(args.output, 'virus.contig.txt')
    # Extract only the ContigID column and write without a header
    contig_df = merged_df[['ContigID']]
    contig_df.to_csv(contig_output_path, sep='\t', index=False, header=False)


if __name__ == "__main__":
    main()