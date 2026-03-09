import argparse
import os
import time
import pandas as pd
from datetime import datetime
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import xml.etree.ElementTree as ET
import urllib.error
import http.client

# simple script that pull all sequences matching a query from genbank. 
# used in background_tree.sh

Entrez.email = "your.email@site.com"

def fetch_reference_sequence(query):
    """Fetch RefSeq reference sequence for a given taxonomic query."""
    search_term = f"{query}[Orgn] AND srcdb_refseq[PROP]"
    try:
        handle = Entrez.esearch(db="nucleotide", term=search_term, sort="relevance", retmax=1)
        record = Entrez.read(handle)
        handle.close()
    except urllib.error.HTTPError:
        time.sleep(5)
        handle = Entrez.esearch(db="nucleotide", term=search_term, sort="relevance", retmax=1)
        record = Entrez.read(handle)
        handle.close()

    if not record["IdList"]:
        raise ValueError(f"No RefSeq record found for {query}.")

    ref_id = record["IdList"][0]
    handle = Entrez.efetch(db="nucleotide", id=ref_id, rettype="gb", retmode="xml")
    xml_string = handle.read()
    handle.close()

    root = ET.fromstring(xml_string)
    gbseq = root.find('GBSeq')
    accession = gbseq.findtext('GBSeq_accession-version')
    length = int(gbseq.findtext('GBSeq_length'))
    sequence = gbseq.findtext('GBSeq_sequence')

    print(f"Fetched reference: {accession} ({length} bp)")
    return accession, length, sequence

def parse_collection_date(raw_date):
    for fmt in ("%d-%b-%Y", "%Y-%m-%d", "%Y", "%b-%Y"):
        try:
            return datetime.strptime(raw_date, fmt).strftime("%Y-%m-%d")
        except ValueError:
            continue
    return ""

def is_human_host(host_value):
    """Return True if the host is acceptable under the --human_only flag."""
    if not host_value:  # empty host field allowed
        return True
    host_lower = host_value.lower()
    return ("human" in host_lower) or ("homo sapiens" in host_lower)

def fetch_and_filter_sequences(search_term, ref_length, ref_accession, min_len_frac, min_gatc_frac, human_only=False):
    Entrez.email = "your.email@site.com"

    # Get total sequence count
    try:
        handle = Entrez.esearch(db="nucleotide", term=search_term, idtype="acc", retmax=0)
        record = Entrez.read(handle)
        handle.close()
        total_count = int(record["Count"])
        print(f"Total sequences found: {total_count}")
    except urllib.error.HTTPError:
        time.sleep(10)
        handle = Entrez.esearch(db="nucleotide", term=search_term, idtype="acc", retmax=0)
        record = Entrez.read(handle)
        handle.close()
        total_count = int(record["Count"])
        print(f"Total sequences found (after retry): {total_count}")

    # Fetch all IDs in batches
    batch_size = 10000
    all_ids = []
    for start in range(0, total_count, batch_size):
        try:
            handle = Entrez.esearch(
                db="nucleotide",
                term=search_term,
                idtype="acc",
                retmax=batch_size,
                retstart=start,
                sort='recently_added'
            )
            record = Entrez.read(handle)
            handle.close()
            all_ids.extend(record["IdList"])
            print(f"Fetched IDs {start} - {start + len(record['IdList'])}")
        except urllib.error.HTTPError:
            time.sleep(10)
            handle = Entrez.esearch(
                db="nucleotide",
                term=search_term,
                idtype="acc",
                retmax=batch_size,
                retstart=start,
                sort='recently_added'
            )
            record = Entrez.read(handle)
            handle.close()
            all_ids.extend(record["IdList"])
            print(f"Fetched IDs {start} - {start + len(record['IdList'])} (after retry)")

    print(f"Total IDs fetched: {len(all_ids)}")

    # Fetch sequences in batches
    seq_records = []
    metadata_dict = {}
    fetch_batch_size = 500  # efetch limit
    for i in range(0, len(all_ids), fetch_batch_size):
        batch_ids = all_ids[i:i+fetch_batch_size]
        try:
            handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="gb", retmode="xml")
            data = handle.read()
            handle.close()
        except urllib.error.HTTPError:
            time.sleep(10)
            handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="gb", retmode="xml")
            data = handle.read()
            handle.close()

        root = ET.fromstring(data)
        for entry in root:
            acc = entry.findtext('.//GBSeq_accession-version')
            seq_text = entry.findtext('.//GBSeq_sequence')
            length = entry.findtext('.//GBSeq_length')
            moltype = entry.findtext('.//GBSeq_moltype')
            primary = entry.findtext('.//GBSeq_primary-accession')

            if not acc or not seq_text:
                continue

            raw_seq = seq_text.upper()
            seq_len = len(raw_seq)
            gatc_fraction = sum(raw_seq.count(b) for b in "GATC") / seq_len if seq_len else 0

            if gatc_fraction < min_gatc_frac or seq_len < min_len_frac * ref_length:
                continue

            meta = {}
            for feat in entry.findall('.//GBSeq_feature-table/GBFeature'):
                for qualifier in feat.findall('./GBFeature_quals/GBQualifier'):
                    key = qualifier.findtext('GBQualifier_name')
                    val = qualifier.findtext('GBQualifier_value')
                    if key and val:
                        meta[key.replace(" ", "_")] = val

            # Optional additions
            meta["accession"] = acc
            meta["primary_accession"] = primary or acc
            meta["length"] = length
            meta["molecule_type"] = moltype or "unknown"

            # Normalize and build description
            date_raw = meta.get("collection_date", "")
            date_iso = parse_collection_date(date_raw)
            meta["collection_date_iso"] = date_iso
            organism = meta.get("organism", "")
            geo = meta.get("geo_loc_name", "")
            host = meta.get("host", "")

            # Apply human filter if requested
            if human_only and not is_human_host(host):
                continue

            desc = " | ".join(filter(None, [organism, geo, host, date_iso]))

            record = SeqRecord(seq=Seq(raw_seq), id=acc, description=desc)
            seq_records.append(record)
            metadata_dict[acc] = meta

    # Sort to place reference first
    seq_records.sort(key=lambda r: (r.id != ref_accession, r.id))
    return seq_records, metadata_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download and filter sequences from NCBI.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pathogen", help="Pathogen name (e.g. 'enterovirus D68')")
    group.add_argument("--taxid", help="NCBI Taxonomy ID (e.g. 42789)")
    parser.add_argument("--min_length_fraction", type=float, default=0.9)
    parser.add_argument("--min_gatc_fraction", type=float, default=0.9)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--human_only", action="store_true", help="Keep only human sequences (host = Human, Homo sapiens, or empty)")

    args = parser.parse_args()

    search_term = f"{args.pathogen}[ORGN]" if args.pathogen else f"txid{args.taxid}[Organism:exp]"
    query_term = args.pathogen or f"txid{args.taxid}"

    ref_acc, ref_len, _ = fetch_reference_sequence(query_term)
    seqs, metadata = fetch_and_filter_sequences(
        search_term,
        ref_len,
        ref_acc,
        args.min_length_fraction,
        args.min_gatc_fraction,
        human_only=args.human_only
    )

    out_dir = os.path.dirname(args.output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    SeqIO.write(seqs, f"{args.output_prefix}.fasta", "fasta")
    pd.DataFrame(metadata).T.to_csv(f"{args.output_prefix}_metadata.csv")
