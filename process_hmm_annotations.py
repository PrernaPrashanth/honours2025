import pandas as pd
import re
import os
import sys

def process_hmm_annotations(df, sample):
    if isinstance(df, str):  # Check if df is a file path
        df = pd.read_csv(df, sep=' ')
    
    df = df.iloc[1: , :]  # Remove the first row
    transcript_dict = {}

    for scaffold, group in df.groupby('#'):
        group = group.copy()  # Avoid SettingWithCopyWarning
        group['subdomain'] = group['Generated'].str.split('_clustalo', expand=True)[0]
        group['start'] = group['one'].astype(str).str.split('-', expand=True)[0].astype(int)

        size_match = re.search(r'size(\d+)', scaffold)
        if not size_match:
            continue  # Skip if no size is found
        size = int(size_match.group(1))

        if size > 1500:
            # Process CIDRa
            count = len(re.findall(r'CIDRa', ','.join(group['subdomain'])))
            if count > 1:
                subset_CIDRa = group[group['Generated'].str.contains('CIDRa')]
                subset_CIDRa['of'] = subset_CIDRa['of'].astype(float)
                only_CIDRa = subset_CIDRa.loc[subset_CIDRa['of'].idxmin()]
                group = group[~group['Generated'].str.contains('CIDRa')]
                group = pd.concat([group, only_CIDRa.to_frame().T]).sort_values(by='start')

            # Process DBLa
            count = len(re.findall(r'DBLa', ','.join(group['subdomain'])))
            if count > 1:
                subset_DBLa = group[group['Generated'].str.contains('DBLa')]
                subset_DBLa['of'] = subset_DBLa['of'].astype(float)
                only_DBLa = subset_DBLa.loc[subset_DBLa['of'].idxmin()]
                group = group[~group['Generated'].str.contains('DBLa')]
                group = pd.concat([group, only_DBLa.to_frame().T]).sort_values(by='start')

            # Construct transcript
            transcript = list(group['subdomain'])
            if len(set(transcript)) == 1:
                transcript = set(transcript)
            transcript = '-'.join(transcript)

            if len(transcript) >= 3:
                transcript_dict[scaffold] = transcript

    # Convert dictionary to DataFrame
    sig_id_ann = pd.DataFrame({'Assembled_id': transcript_dict.keys(), 'Annotation': transcript_dict.values()})
    file_name = f"{sample}_assembled_id_and_var_annotation.csv"
    sig_id_ann.to_csv(file_name, index=False)

    # Save individual files for each scaffold with domain info
    for scaffold, transcript in transcript_dict.items():
        # Extract domain type from the transcript
        domain_type = transcript.replace('-', '_')  # Convert hyphens to underscores for filenames
        
        # Create base filename
        base_filename = f"{sample}_{scaffold}_{domain_type}"
        filename = f"{base_filename}.txt"
        
        # Check if file already exists and add suffix if needed
        suffix = 1
        while os.path.exists(filename):
            filename = f"{base_filename}_{suffix}.txt"
            suffix += 1
        
        # Create DataFrame with just this ID and save
        single_id_df = pd.DataFrame({'Assembled_id': [scaffold]})
        single_id_df.to_csv(filename, index=False)

    # Still save the original combined file
    sig_id_only = pd.DataFrame({'Assembled_id': sig_id_ann['Assembled_id']})
    file_name = f"{sample}_assembled_id_sig_annotation.txt" 
    sig_id_only.to_csv(file_name, index=False)

    if len(transcript_dict) == 0:
        print(f"Note: No transcripts found for {sample} — all scaffolds likely <=1500 bp or filtered out.")

    return transcript_dict

# Add code to run the function when script is executed
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_hmm.py input_file sample_name")
        sys.exit(1)
    
    input_file = sys.argv[1]
    sample_name = sys.argv[2]
    
    result = process_hmm_annotations(input_file, sample_name)
    print(f"Processing complete. Found {len(result)} transcripts.")


