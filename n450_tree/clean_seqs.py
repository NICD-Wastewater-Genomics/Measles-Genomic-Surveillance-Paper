#!/usr/bin/env python3
import sys

## removes individual sequences and shared regions with all Ns (e.g. trim to only N450)
def read_fasta(filename):
    headers = []
    seqs = []
    seq = []

    with open(filename) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                headers.append(line)
                if seq:
                    seqs.append("".join(seq))
                    seq = []
            else:
                seq.append(line)
        if seq:
            seqs.append("".join(seq))

    return headers, seqs


def write_fasta(headers, seqs, out):
    with open(out, "w") as fh:
        for h, s in zip(headers, seqs):
            fh.write(h + "\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i+80] + "\n")


def remove_all_N_columns(headers, seqs):
    nseq = len(seqs)
    L = len(seqs[0])

    # Identify columns to keep (at least one non-N across seqs)
    keep = []
    for i in range(L):
        col = [seqs[s][i] for s in range(nseq)]
        if not all(c.upper() == "N" for c in col):
            keep.append(i)

    # Filter columns
    filtered = []
    for s in seqs:
        filtered.append("".join(s[i] for i in keep))

    return filtered


def remove_all_N_sequences(headers, seqs):
    """Remove sequences that are all Ns after column filtering."""
    new_headers = []
    new_seqs = []
    for h, s in zip(headers, seqs):
        if not all(c.upper() == "N" for c in s):
            new_headers.append(h)
            new_seqs.append(s)
    return new_headers, new_seqs


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(1)

    inp, out = sys.argv[1], sys.argv[2]

    headers, seqs = read_fasta(inp)

    # Sanity: aligned?
    lengths = set(len(s) for s in seqs)
    if len(lengths) != 1:
        print("Error: sequences have different lengths (not aligned).")
        sys.exit(1)

    # Remove all-N columns
    seqs = remove_all_N_columns(headers, seqs)

    # Remove sequences that are now all Ns
    headers, seqs = remove_all_N_sequences(headers, seqs)

    # Write result
    write_fasta(headers, seqs, out)

