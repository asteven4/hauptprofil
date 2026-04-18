import pandas as pd
import re

# Load the clinical dataset
clinical_df = pd.read_csv(
    "C:/Users/astev/PycharmProjects/tsvmerger/selected_dataset_mk14.tsv",
    sep="\t",
    encoding="latin1",
    low_memory=False
)

# Extract patient barcodes from the clinical dataset
patient_barcodes = set(clinical_df["bcr_patient_barcode"])

# Large gene expression dataset
input_file = "C:/Users/astev/PycharmProjects/tsvmerger/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"
output_file = "C:/Users/astev/PycharmProjects/tsvmerger_mk2/reformatted_gene_df.tsv"

# Read and process the large gene file in chunks to avoid memory overload
chunk_size = 30000
reader = pd.read_csv(input_file, sep='\t', index_col=0, chunksize=chunk_size)

# Open the output file once to write chunks iteratively
with open(output_file, 'w', encoding="utf-8") as f_out:
    for i, chunk in enumerate(reader):
        chunk = chunk.T.reset_index()  # Transpose, make index a column
        chunk.rename(columns={"index": "bcr_patient_barcode"}, inplace=True)

        # Append mode writing to handle large files efficiently
        chunk.to_csv(f_out, sep='\t', index=False, header=(i == 0))

print("Large file processed in chunks and saved as:", output_file)

# Reload the processed file with corrected formatting
file_to_be_cleaned = pd.read_csv(output_file, sep="\t", low_memory=False)


# Function to clean patient IDs by removing extra ID parts
def clean_column_name(col_name):
    match = re.match(r"([A-Z0-9-]+)", col_name)  # Extracts the main barcode before extra ID parts
    return match.group(1) if match else col_name  # Return cleaned ID


# Apply transformation to clean the column names
file_to_be_cleaned["bcr_patient_barcode"] = file_to_be_cleaned["bcr_patient_barcode"].str.extract(r"([A-Z0-9-]{12})")[0]

# Filter the columns to only those matching patient barcodes
# Keep all gene expression columns while selecting only relevant patients
filtered_expression_df = file_to_be_cleaned[file_to_be_cleaned["bcr_patient_barcode"].isin(patient_barcodes)]
# Merge the cleaned gene expression data with clinical data
merged_df = clinical_df.merge(filtered_expression_df, on="bcr_patient_barcode", how="inner")

# Save the final merged dataset
merged_output_file = "C:/Users/astev/PycharmProjects/tsvmerger_mk2/MergedClinical+PANCAN_mk3.tsv"
merged_df.to_csv(merged_output_file, sep="\t", index=False)

print("Processing complete. Merged file saved as:", merged_output_file)
