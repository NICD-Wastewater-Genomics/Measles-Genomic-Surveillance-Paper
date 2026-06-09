#!/usr/bin/env python3

from Bio import SeqIO
import csv
import pandas as pd

MSA_FILE = "aligned.fasta"
BED_FILE = "primer_v3_400.bed"
OUTPUT_CSV = "primer_sequences.csv"

# 1. Read and parse the BED file using the pandas approach
col_names = [
    "ref",
    "start",
    "end",
    "left_right",
    "primer_pool",
    "strand",
    "primer_seq",
]

primer_bed = pd.read_csv(BED_FILE, sep="\t", names=col_names, comment="#")

# Parse out the core amplicon number from the name (e.g., "Primer_1_LEFT" -> "1")
primer_bed["amplicon_number"] = primer_bed["left_right"].str.split("_").str[1]

# 2. Merge the DataFrame with itself to pair LEFT and RIGHT primers on one row
df_paired = pd.merge(
    primer_bed.loc[primer_bed["left_right"].str.contains("LEFT")],
    primer_bed.loc[primer_bed["left_right"].str.contains("RIGHT")],
    on=["amplicon_number", "primer_pool"],
    suffixes=("_left", "_right"),
)

# Handle duplicate amplicon numbers if they exist (e.g., alternative/overlap primers)
mask = df_paired.duplicated(subset=["amplicon_number"], keep=False)
if mask.any():
    df_paired.loc[mask, "amplicon_number"] = (
        df_paired.loc[mask, "amplicon_number"].astype(str)
        + "_"
        + df_paired.loc[mask].groupby("amplicon_number").cumcount().astype(str)
    )

# 3. Process the MSA and extract sequence segments for each amplicon pair
# Load MSA records into memory
records = list(SeqIO.parse(MSA_FILE, "fasta"))

output_rows = []

for rec in records:
    seq = str(rec.seq)

    # Iterate through the paired amplicons
    for _, row in df_paired.iterrows():
        # Extract sequences from the MSA using the coordinates
        left_extracted = seq[int(row["start_left"]) : int(row["end_left"])]
        right_extracted = seq[int(row["start_right"]) : int(row["end_right"])]

        output_rows.append(
            {
                "sequence_id": rec.id,
                "amplicon_number": row["amplicon_number"],
                "pool": row["primer_pool"],
                "left_primer_name": row["left_right_left"],
                "left_primer_seq_ref": row["primer_seq_left"],
                "left_aligned_sequence": left_extracted,
                "right_primer_name": row["left_right_right"],
                "right_primer_seq_ref": row["primer_seq_right"],
                "right_aligned_sequence": right_extracted,
            }
        )

# 4. Save results to a clean CSV
output_df = pd.DataFrame(output_rows)
output_df.to_csv(OUTPUT_CSV, index=False)

print(f"Output successfully written to {OUTPUT_CSV}")