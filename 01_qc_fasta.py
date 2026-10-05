#!/usr/bin/env python3

import argparse
import os
import sys
from collections import Counter

from Bio import SeqIO


# Standard 20 amino acids plus common ambiguous/unknown amino-acid codes.
ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Perform basic quality control on a protein FASTA file "
            "before structural phylogenetic analysis."
        )
    )

    parser.add_argument(
        "--input_fasta",
        required=True,
        help="Full path to the input protein FASTA file."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help=(
            "Full path to the output directory. "
            "A qc/ subdirectory will be created here."
        )
    )

    return parser.parse_args()


def main():

    args = parse_args()

    input_fasta = os.path.abspath(args.input_fasta)
    output_dir = os.path.abspath(args.output_dir)

    # ---------------------------------------------------------
    # Check input
    # ---------------------------------------------------------

    if not os.path.isfile(input_fasta):
        sys.exit(f"ERROR: Input FASTA does not exist:\n{input_fasta}")

    # ---------------------------------------------------------
    # Set up output
    # ---------------------------------------------------------

    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)

    output_tsv = os.path.join(qc_dir, "sequence_qc.tsv")

    # ---------------------------------------------------------
    # Read FASTA
    # ---------------------------------------------------------

    records = list(SeqIO.parse(input_fasta, "fasta"))

    if len(records) == 0:
        sys.exit(f"ERROR: No FASTA records found in:\n{input_fasta}")

    # ---------------------------------------------------------
    # Check duplicate IDs
    # ---------------------------------------------------------

    ids = [record.id for record in records]
    id_counts = Counter(ids)

    duplicate_ids = {
        seq_id
        for seq_id, count in id_counts.items()
        if count > 1
    }

    # ---------------------------------------------------------
    # Analyze sequences
    # ---------------------------------------------------------

    qc_rows = []

    for record in records:

        seq_id = record.id
        sequence = str(record.seq).upper()

        invalid_chars = sorted(set(sequence) - ALLOWED_AA)

        duplicate = seq_id in duplicate_ids
        valid_characters = len(invalid_chars) == 0

        qc_rows.append({
            "sequence_id": seq_id,
            "length": len(sequence),
            "valid_characters": valid_characters,
            "invalid_characters": ",".join(invalid_chars),
            "duplicate_id": duplicate
        })

    # ---------------------------------------------------------
    # Write QC table
    # ---------------------------------------------------------

    with open(output_tsv, "w") as outfile:

        outfile.write(
            "sequence_id\t"
            "length\t"
            "valid_characters\t"
            "invalid_characters\t"
            "duplicate_id\n"
        )

        for row in qc_rows:

            outfile.write(
                f"{row['sequence_id']}\t"
                f"{row['length']}\t"
                f"{row['valid_characters']}\t"
                f"{row['invalid_characters']}\t"
                f"{row['duplicate_id']}\n"
            )

    # ---------------------------------------------------------
    # Summary statistics
    # ---------------------------------------------------------

    lengths = [row["length"] for row in qc_rows]

    invalid_sequences = [
        row for row in qc_rows
        if not row["valid_characters"]
    ]

    print()
    print("========================================")
    print("Protein FASTA QC")
    print("========================================")
    print(f"Input FASTA:       {input_fasta}")
    print(f"Output directory:  {output_dir}")
    print()
    print(f"Sequences:         {len(records)}")
    print(f"Shortest:          {min(lengths)} aa")
    print(f"Longest:           {max(lengths)} aa")
    print(f"Mean length:       {sum(lengths) / len(lengths):.1f} aa")
    print()

    if duplicate_ids:
        print(f"WARNING: {len(duplicate_ids)} duplicate sequence IDs found:")
        for seq_id in sorted(duplicate_ids):
            print(f"  {seq_id}")
    else:
        print("Duplicate IDs:     None")

    if invalid_sequences:
        print(
            f"WARNING: {len(invalid_sequences)} sequences "
            "contain invalid characters."
        )

        for row in invalid_sequences:
            print(
                f"  {row['sequence_id']}: "
                f"{row['invalid_characters']}"
            )
    else:
        print("Invalid sequences: None")

    print()
    print(f"QC table written:  {output_tsv}")
    print("========================================")
    print()


if __name__ == "__main__":
    main()