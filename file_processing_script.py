import pandas as pd
import numpy as np
import glob
import os

# Read the Excel sheet with disease classification
# Replace 'your_excel_file.xlsx' with your actual Excel file name
excel_df = pd.read_excel('Book1.xlsx')  # Adjust sheet name if needed

# Find all files matching the pattern *_transcripts_count.txt
file_pattern = "*_transcript_counts.txt"
transcript_files = glob.glob(file_pattern)

# Extract sample names by removing the "_transcripts_count.txt" suffix
files = [os.path.basename(f).replace('_transcript_counts.txt', '') for f in transcript_files]

print(f"Found {len(files)} transcript count files")

# Separate CH and non-CH files
ch_files = [f for f in files if f.startswith('CH')]
non_ch_files = [f for f in files if not f.startswith('CH')]

# Create a mapping function to convert CH file names to Excel format
def ch_to_excel_format(ch_name):
    """Convert CH1785 format to 1785-CH format"""
    number = ch_name[2:]  # Remove 'CH' prefix
    return f"{number}-CH"

# Create a dictionary mapping Excel format to original CH file names
ch_mapping = {ch_to_excel_format(ch): ch for ch in ch_files}

# Get the disease classification column from Excel
# Using the correct column names from your Excel file
sample_col = 'Sample_Nr'  # Column with sample names in Excel (format: 1785-CH)
disease_col = 'DISEASE'  # Column with disease classification

# Filter Excel data for CH samples only
ch_excel_data = excel_df[excel_df[sample_col].str.contains('-CH', na=False)].copy()

# Sort CH samples by disease classification
ch_excel_data_sorted = ch_excel_data.sort_values(by=disease_col)

# Create ordered list of CH files based on disease classification
ordered_ch_files = []
for excel_sample in ch_excel_data_sorted[sample_col]:
    if excel_sample in ch_mapping:
        ordered_ch_files.append(ch_mapping[excel_sample])

# Add any CH files not found in Excel to the end
missing_ch_files = [ch for ch in ch_files if ch not in ordered_ch_files]
ordered_ch_files.extend(missing_ch_files)

# Final ordered file list: non-CH files first, then sorted CH files
final_file_order = non_ch_files + ordered_ch_files

print("File processing order:")
for i, file in enumerate(final_file_order, 1):
    print(f"{i}. {file}")

# Process files in the new order
df_list = []
for file in final_file_order:
    colnames = ('geneid', file+'_counts')
    df = pd.read_csv(file+'_transcript_counts.txt', sep='\t', engine='python', names=colnames)
    df_list.append(df)

# Merge dataframes
merged_df = df_list[0]
for i in range(1, len(df_list)):
    merged_df = pd.merge(merged_df, df_list[i], on='geneid', how='outer')

# Fill NaN values with 0 and round
merged_df = merged_df.fillna(0)
merged_df = merged_df.round(0)

# Convert to int to remove decimal points
for col in merged_df.columns:
    if col != 'geneid':
        merged_df[col] = merged_df[col].astype(int)

# Save the result
merged_df.to_csv("SM_domain_counts.txt", sep='\t', index=False)

print(f"\nProcessing complete! Output saved to SM_domain_counts.txt")
print(f"Final dataframe shape: {merged_df.shape}")
print(f"Columns: {list(merged_df.columns)}")

# Create a mapping of CH file name -> disease classification
ch_disease_map = {}
for _, row in ch_excel_data.iterrows():
    excel_name = row[sample_col]      # e.g. "1785-CH"
    disease = row[disease_col]
    if excel_name in ch_mapping:
        ch_disease_map[ch_mapping[excel_name]] = disease

# Build a list of tuples (sample_name, disease_classification_or_blank)
sample_order_with_disease = []
for sample in final_file_order:
    disease = ch_disease_map.get(sample, "")  # Empty if not CH or missing
    sample_order_with_disease.append((sample, disease))

# Show the table
print("\nSample order with disease classification:")
for i, (sample, disease) in enumerate(sample_order_with_disease, 1):
    print(f"{i}. {sample}\t{disease}")

# Optional: Save to file
pd.DataFrame(sample_order_with_disease, columns=["Sample", "Disease"]).to_csv(
    "sample_order_with_disease.txt", sep="\t", index=False
)

