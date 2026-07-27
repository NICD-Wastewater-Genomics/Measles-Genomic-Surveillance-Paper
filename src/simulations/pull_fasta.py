import pandas as pd
from Bio import Entrez, SeqIO
import argparse

def fetch_fasta(accessions, output_file):
    Entrez.email = "X"

    with Entrez.efetch(
        db="nucleotide",
        id=accessions,
        rettype="fasta",
        retmode="text"
    ) as handle:

        records = list(SeqIO.parse(handle, "fasta"))

        # Deduplicate by the first word of the FASTA header
        seen = set()
        unique_records = []

        for record in records:
            seq_id = record.description.split()[0]

            if seq_id not in seen:
                seen.add(seq_id)
                record.description = ""  # optional
                unique_records.append(record)

        with open(output_file, "w") as out_handle:
            SeqIO.write(unique_records, out_handle, "fasta")


def main():
    parser = argparse.ArgumentParser(
        description="Pull fasta files from a list of accessions"
    )
    parser.add_argument(
        "accession_list",
        type=str,
        help="Input TSV file containing accession list"
    )
    parser.add_argument(
        "output_fasta",
        type=str,
        help="Output FASTA file for the downloaded sequences"
    )

    args = parser.parse_args()

    metadata_df = pd.read_csv(args.accession_list, sep="\t")
    accessions = metadata_df["V1"].tolist()

    fetch_fasta(accessions, args.output_fasta)


if __name__ == "__main__":
    main()