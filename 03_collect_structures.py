#!/usr/bin/env python3

import argparse
import glob
import json
import os
import shutil
import sys

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Collect rank-1 ColabFold structures and confidence scores "
            "for structural phylogenetic analysis."
        )
    )

    parser.add_argument(
        "--input_fasta",
        required=True,
        help="Full path to the original protein FASTA file."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help=(
            "Full path to the analysis output directory containing "
            "colabfold_results/."
        )
    )

    return parser.parse_args()


def find_single_file(pattern, description, sequence_id):
    """
    Find exactly one file matching a glob pattern.
    Return the path if exactly one exists.
    Otherwise report an error.
    """

    matches = glob.glob(pattern)

    if len(matches) == 0:
        print(
            f"ERROR: No {description} found for {sequence_id}",
            file=sys.stderr
        )
        return None

    if len(matches) > 1:
        print(
            f"ERROR: Multiple {description} files found for {sequence_id}:",
            file=sys.stderr
        )

        for match in matches:
            print(f"  {match}", file=sys.stderr)

        return None

    return matches[0]


def main():

    args = parse_args()

    input_fasta = os.path.abspath(args.input_fasta)
    output_dir = os.path.abspath(args.output_dir)

    colabfold_dir = os.path.join(output_dir, "colabfold_results")
    structure_dir = os.path.join(output_dir, "structures")
    qc_dir = os.path.join(output_dir, "qc")

    output_tsv = os.path.join(qc_dir, "structure_qc.tsv")

    # ---------------------------------------------------------
    # Check inputs
    # ---------------------------------------------------------

    if not os.path.isfile(input_fasta):
        sys.exit(
            f"ERROR: Input FASTA does not exist:\n{input_fasta}"
        )

    if not os.path.isdir(colabfold_dir):
        sys.exit(
            "ERROR: ColabFold results directory does not exist:\n"
            f"{colabfold_dir}"
        )

    # ---------------------------------------------------------
    # Read FASTA
    # ---------------------------------------------------------

    records = list(SeqIO.parse(input_fasta, "fasta"))

    if len(records) == 0:
        sys.exit(
            f"ERROR: No FASTA records found in:\n{input_fasta}"
        )

    # ---------------------------------------------------------
    # IMPORTANT:
    # First validate the COMPLETE ColabFold run.
    #
    # We do this before creating/copying structures so that a
    # partial ColabFold run cannot produce a misleading partial
    # structures/ directory.
    # ---------------------------------------------------------

    collected = []
    errors = 0

    print()
    print("========================================")
    print("Validating ColabFold results")
    print("========================================")
    print(f"Input FASTA:       {input_fasta}")
    print(f"ColabFold results: {colabfold_dir}")
    print(f"Expected proteins: {len(records)}")
    print()

    for record in records:

        sequence_id = record.id
        sequence_length = len(record.seq)

        # Check ColabFold completion marker

        done_file = os.path.join(
            colabfold_dir,
            f"{sequence_id}.done.txt"
        )

        if not os.path.isfile(done_file):
            print(
                f"ERROR: Missing completion marker for {sequence_id}"
            )
            errors += 1
            continue

        # Find rank-1 PDB

        pdb_pattern = os.path.join(
            colabfold_dir,
            f"{sequence_id}_unrelaxed_rank_001_*.pdb"
        )

        pdb_file = find_single_file(
            pdb_pattern,
            "rank-1 PDB",
            sequence_id
        )

        if pdb_file is None:
            errors += 1
            continue

        # Find rank-1 score JSON

        json_pattern = os.path.join(
            colabfold_dir,
            f"{sequence_id}_scores_rank_001_*.json"
        )

        json_file = find_single_file(
            json_pattern,
            "rank-1 score JSON",
            sequence_id
        )

        if json_file is None:
            errors += 1
            continue

        # Read confidence scores

        try:
            with open(json_file, "r") as infile:
                scores = json.load(infile)
        except Exception as exc:
            print(
                f"ERROR: Could not read JSON for {sequence_id}: {exc}"
            )
            errors += 1
            continue

        # ColabFold rank score JSON contains per-residue pLDDT
        # values and a scalar pTM value.

        plddt_values = scores.get("plddt")
        ptm = scores.get("ptm")

        if plddt_values is None:
            print(
                f"ERROR: No pLDDT values found for {sequence_id}"
            )
            errors += 1
            continue

        if ptm is None:
            print(
                f"ERROR: No pTM value found for {sequence_id}"
            )
            errors += 1
            continue

        mean_plddt = sum(plddt_values) / len(plddt_values)

        high_confidence = sum(
            value >= 90 for value in plddt_values
        )

        low_confidence = sum(
            value < 50 for value in plddt_values
        )

        percent_plddt_ge_90 = (
            100 * high_confidence / len(plddt_values)
        )

        percent_plddt_lt_50 = (
            100 * low_confidence / len(plddt_values)
        )

        collected.append({
            "sequence_id": sequence_id,
            "length": sequence_length,
            "mean_plddt": mean_plddt,
            "ptm": ptm,
            "percent_plddt_ge_90": percent_plddt_ge_90,
            "percent_plddt_lt_50": percent_plddt_lt_50,
            "pdb_file": pdb_file,
            "json_file": json_file
        })

    # ---------------------------------------------------------
    # Stop if anything is missing
    # ---------------------------------------------------------

    if errors > 0 or len(collected) != len(records):

        print()
        print("========================================")
        print("VALIDATION FAILED")
        print("========================================")
        print(f"Expected proteins: {len(records)}")
        print(f"Valid predictions: {len(collected)}")
        print(f"Errors:            {errors}")
        print()
        print(
            "No structures were copied. "
            "Check the ColabFold results above."
        )
        print("========================================")
        print()

        sys.exit(1)

    # ---------------------------------------------------------
    # Create output directories only after validation succeeds
    # ---------------------------------------------------------

    os.makedirs(structure_dir, exist_ok=True)
    os.makedirs(qc_dir, exist_ok=True)

    # ---------------------------------------------------------
    # Copy and rename rank-1 structures
    # ---------------------------------------------------------

    for row in collected:

        destination = os.path.join(
            structure_dir,
            f"{row['sequence_id']}.pdb"
        )

        shutil.copy2(
            row["pdb_file"],
            destination
        )

    # ---------------------------------------------------------
    # Write structure QC table
    # ---------------------------------------------------------

    with open(output_tsv, "w") as outfile:

        outfile.write(
            "sequence_id\t"
            "length\t"
            "mean_plddt\t"
            "ptm\t"
            "percent_plddt_ge_90\t"
            "percent_plddt_lt_50\n"
        )

        for row in collected:

            outfile.write(
                f"{row['sequence_id']}\t"
                f"{row['length']}\t"
                f"{row['mean_plddt']:.2f}\t"
                f"{row['ptm']:.4f}\t"
                f"{row['percent_plddt_ge_90']:.2f}\t"
                f"{row['percent_plddt_lt_50']:.2f}\n"
            )

    # ---------------------------------------------------------
    # Summary
    # ---------------------------------------------------------

    mean_family_plddt = (
        sum(row["mean_plddt"] for row in collected)
        / len(collected)
    )

    mean_family_ptm = (
        sum(row["ptm"] for row in collected)
        / len(collected)
    )

    print()
    print("========================================")
    print("Structure collection complete")
    print("========================================")
    print(f"Proteins:          {len(collected)}")
    print(f"Structures copied: {len(collected)}")
    print()
    print(f"Mean family pLDDT: {mean_family_plddt:.2f}")
    print(f"Mean family pTM:   {mean_family_ptm:.3f}")
    print()
    print(f"Structures:        {structure_dir}")
    print(f"QC table:          {output_tsv}")
    print("========================================")
    print()


if __name__ == "__main__":
    main()