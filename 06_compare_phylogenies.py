#!/usr/bin/env python3

"""
06_compare_phylogenies.py

Compare standard amino-acid and structural (3Di) phylogenies.

This script:

1. Reads the standard and structural IQ-TREE trees with Biopython.
2. Reads a TSV containing species prefixes representing a convergent trait.
3. Identifies convergent taxa in each tree based on species prefix.
4. Plots each tree in an unrooted layout:
       - convergent taxa shown in red
       - other taxa shown in black
       - bootstrap support shown at internal nodes
5. Calculates all pairwise topology-only distances among convergent tips.
       - distance = number of edges separating two tips
       - branch lengths are ignored
6. Writes pairwise-distance tables.
7. Produces a KDE plot comparing the pairwise-distance distributions.

Usage
-----

python 06_compare_phylogenies.py \
    --standard_tree /path/to/standard.treefile \
    --structural_tree /path/to/structural.treefile \
    --trait_tsv /path/to/convergent_trait.tsv \
    --output_dir /path/to/output_directory
"""

import argparse
import csv
import math
import os
from collections import deque
from itertools import combinations

import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo


# ============================================================
# ARGUMENTS
# ============================================================

parser = argparse.ArgumentParser(
    description="Plot and compare standard and structural phylogenies."
)

parser.add_argument(
    "--standard_tree",
    required=True,
    help="Full path to standard amino-acid IQ-TREE treefile."
)

parser.add_argument(
    "--structural_tree",
    required=True,
    help="Full path to structural 3Di IQ-TREE treefile."
)

parser.add_argument(
    "--trait_tsv",
    required=True,
    help="TSV containing column 'Convergent_species'."
)

parser.add_argument(
    "--output_dir",
    required=True,
    help="Directory for plots and distance tables."
)

args = parser.parse_args()


# ============================================================
# PATHS
# ============================================================

standard_tree_path = os.path.abspath(args.standard_tree)
structural_tree_path = os.path.abspath(args.structural_tree)
trait_tsv_path = os.path.abspath(args.trait_tsv)
output_dir = os.path.abspath(args.output_dir)

os.makedirs(output_dir, exist_ok=True)


# ============================================================
# INPUT VALIDATION
# ============================================================

for path in [
    standard_tree_path,
    structural_tree_path,
    trait_tsv_path
]:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")


# ============================================================
# READ TREES
# ============================================================

print()
print("============================================================")
print("Reading trees")
print("============================================================")
print()

standard_tree = Phylo.read(standard_tree_path, "newick")
structural_tree = Phylo.read(structural_tree_path, "newick")

standard_tips = standard_tree.get_terminals()
structural_tips = structural_tree.get_terminals()

print(f"Standard tree taxa:   {len(standard_tips)}")
print(f"Structural tree taxa: {len(structural_tips)}")


# ============================================================
# VERIFY TAXON SETS
# ============================================================

standard_names = {tip.name for tip in standard_tips}
structural_names = {tip.name for tip in structural_tips}

if standard_names == structural_names:

    print("Taxon sets: IDENTICAL")

else:

    print()
    print("WARNING: Taxon sets are NOT identical.")

    only_standard = sorted(standard_names - structural_names)
    only_structural = sorted(structural_names - standard_names)

    if only_standard:
        print()
        print("Only in standard tree:")
        for name in only_standard:
            print(f"  {name}")

    if only_structural:
        print()
        print("Only in structural tree:")
        for name in only_structural:
            print(f"  {name}")


# ============================================================
# READ CONVERGENT SPECIES PREFIXES
# ============================================================

print()
print("============================================================")
print("Reading convergent-trait species")
print("============================================================")
print()

convergent_prefixes = []

with open(trait_tsv_path, newline="") as handle:

    reader = csv.DictReader(handle, delimiter="\t")

    if reader.fieldnames is None or "Convergent_species" not in reader.fieldnames:
        raise ValueError(
            "Trait TSV must contain a column named 'Convergent_species'."
        )

    for row in reader:

        prefix = row["Convergent_species"].strip()

        if prefix:
            convergent_prefixes.append(prefix)


convergent_prefixes = set(convergent_prefixes)

print(f"Convergent species prefixes: {len(convergent_prefixes)}")

for prefix in sorted(convergent_prefixes):
    print(f"  {prefix}")


# ============================================================
# IDENTIFY CONVERGENT TIPS
# ============================================================

def species_prefix(tip_name):
    """
    Extract species prefix from tree tip.

    Example:
        DroSpa_Ds00010569RA -> DroSpa
    """
    return tip_name.split("_", 1)[0]


def get_convergent_tips(tree):

    tips = []

    for tip in tree.get_terminals():

        prefix = species_prefix(tip.name)

        if prefix in convergent_prefixes:
            tips.append(tip)

    return tips


standard_convergent = get_convergent_tips(standard_tree)
structural_convergent = get_convergent_tips(structural_tree)


print()
print(f"Convergent tips in standard tree:   {len(standard_convergent)}")
print(f"Convergent tips in structural tree: {len(structural_convergent)}")


print()
print("Convergent tips in standard tree:")

for tip in standard_convergent:
    print(f"  {tip.name}")


print()
print("Convergent tips in structural tree:")

for tip in structural_convergent:
    print(f"  {tip.name}")


# ============================================================
# BUILD UNDIRECTED TREE GRAPH
# ============================================================

def build_adjacency(tree):
    """
    Convert a Biopython tree into an undirected adjacency dictionary.

    Branch lengths are deliberately ignored.
    """

    adjacency = {}

    for parent in tree.find_clades(order="level"):

        adjacency.setdefault(parent, [])

        for child in parent.clades:

            adjacency.setdefault(child, [])

            adjacency[parent].append(child)
            adjacency[child].append(parent)

    return adjacency


# ============================================================
# TOPOLOGICAL DISTANCE
# ============================================================

def edge_distance(adjacency, start, end):
    """
    Number of edges on the unique path between two tree nodes.

    This is a topology-only distance.
    Branch lengths are ignored.
    """

    queue = deque([(start, 0)])
    visited = {start}

    while queue:

        node, distance = queue.popleft()

        if node is end:
            return distance

        for neighbor in adjacency[node]:

            if neighbor not in visited:

                visited.add(neighbor)
                queue.append((neighbor, distance + 1))

    raise RuntimeError("No path found between tree nodes.")


# ============================================================
# PAIRWISE CONVERGENT-TIP DISTANCES
# ============================================================

def calculate_pairwise_distances(tree, convergent_tips):

    adjacency = build_adjacency(tree)

    results = []

    for tip1, tip2 in combinations(convergent_tips, 2):

        distance = edge_distance(
            adjacency,
            tip1,
            tip2
        )

        results.append({
            "tip1": tip1.name,
            "tip2": tip2.name,
            "species1": species_prefix(tip1.name),
            "species2": species_prefix(tip2.name),
            "topological_distance": distance
        })

    return results


standard_distances = calculate_pairwise_distances(
    standard_tree,
    standard_convergent
)

structural_distances = calculate_pairwise_distances(
    structural_tree,
    structural_convergent
)


# ============================================================
# WRITE DISTANCE TABLES
# ============================================================

def write_distance_table(results, output_file):

    fieldnames = [
        "tip1",
        "tip2",
        "species1",
        "species2",
        "topological_distance"
    ]

    with open(output_file, "w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(results)


standard_distance_file = os.path.join(
    output_dir,
    "standard_convergent_pairwise_distances.tsv"
)

structural_distance_file = os.path.join(
    output_dir,
    "structural_convergent_pairwise_distances.tsv"
)


write_distance_table(
    standard_distances,
    standard_distance_file
)

write_distance_table(
    structural_distances,
    structural_distance_file
)


# ============================================================
# UNROOTED TREE LAYOUT
# ============================================================

def initial_circle_layout(tree):
    """
    Construct an initial unrooted layout.

    Tips are distributed evenly around a circle.
    Internal nodes begin at the mean position of neighboring
    descendant tips.
    """

    terminals = tree.get_terminals()

    positions = {}

    n = len(terminals)

    # Tips around a circle
    for i, tip in enumerate(terminals):

        angle = 2 * math.pi * i / n

        positions[tip] = np.array([
            math.cos(angle),
            math.sin(angle)
        ], dtype=float)


    # Internal nodes positioned at mean descendant-tip position
    for clade in reversed(list(tree.find_clades(order="level"))):

        if clade.is_terminal():
            continue

        descendants = clade.get_terminals()

        coords = np.array([
            positions[x]
            for x in descendants
        ])

        positions[clade] = coords.mean(axis=0)

    return positions


def relax_layout(tree, positions, iterations=400):
    """
    Relax internal node positions while keeping tips fixed.

    This produces an unrooted topology-oriented visualization.
    Branch lengths are not used.
    """

    adjacency = build_adjacency(tree)

    terminals = set(tree.get_terminals())

    internal_nodes = [
        node
        for node in adjacency
        if node not in terminals
    ]

    for _ in range(iterations):

        new_positions = {}

        for node in internal_nodes:

            neighbors = adjacency[node]

            coords = np.array([
                positions[n]
                for n in neighbors
            ])

            new_positions[node] = coords.mean(axis=0)

        for node, pos in new_positions.items():
            positions[node] = pos

    return positions


# ============================================================
# TREE PLOTTING
# ============================================================

def plot_unrooted_tree(tree, output_file, title):

    adjacency = build_adjacency(tree)

    positions = initial_circle_layout(tree)

    positions = relax_layout(
        tree,
        positions
    )

    fig, ax = plt.subplots(
        figsize=(16, 16)
    )


    # --------------------------------------------------------
    # Draw edges
    # --------------------------------------------------------

    drawn_edges = set()

    for node, neighbors in adjacency.items():

        for neighbor in neighbors:

            edge_key = frozenset([node, neighbor])

            if edge_key in drawn_edges:
                continue

            drawn_edges.add(edge_key)

            x1, y1 = positions[node]
            x2, y2 = positions[neighbor]

            ax.plot(
                [x1, x2],
                [y1, y2],
                linewidth=0.8,
                color="black",
                zorder=1
            )


    # --------------------------------------------------------
    # Draw bootstrap support
    # --------------------------------------------------------

    for clade in tree.get_nonterminals():

        if clade.confidence is None:
            continue

        x, y = positions[clade]

        support = float(clade.confidence)

        if support.is_integer():
            support_text = str(int(support))
        else:
            support_text = f"{support:.1f}"

        ax.text(
            x,
            y,
            support_text,
            fontsize=5.5,
            ha="center",
            va="center",
            bbox=dict(
                facecolor="white",
                edgecolor="none",
                alpha=0.75,
                pad=0.5
            ),
            zorder=3
        )


    # --------------------------------------------------------
    # Draw tip labels
    # --------------------------------------------------------

    for tip in tree.get_terminals():

        x, y = positions[tip]

        prefix = species_prefix(tip.name)

        is_convergent = prefix in convergent_prefixes

        label_color = "red" if is_convergent else "black"

        # Offset label away from center
        vector = np.array([x, y])

        norm = np.linalg.norm(vector)

        if norm > 0:
            vector = vector / norm

        label_x = x + vector[0] * 0.035
        label_y = y + vector[1] * 0.035

        if vector[0] >= 0:
            ha = "left"
        else:
            ha = "right"

        ax.text(
            label_x,
            label_y,
            tip.name,
            fontsize=6.5,
            color=label_color,
            ha=ha,
            va="center",
            zorder=4
        )


    # --------------------------------------------------------
    # Formatting
    # --------------------------------------------------------

    ax.set_title(
        title,
        fontsize=16,
        pad=20
    )

    ax.set_aspect("equal")

    ax.axis("off")

    # Extra room for long taxon labels
    ax.margins(0.30)

    plt.tight_layout()

    plt.savefig(
        output_file,
        format="pdf",
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# CREATE TREE FIGURES
# ============================================================

print()
print("============================================================")
print("Plotting trees")
print("============================================================")
print()


standard_tree_pdf = os.path.join(
    output_dir,
    "standard_tree_unrooted.pdf"
)

structural_tree_pdf = os.path.join(
    output_dir,
    "structural_3di_tree_unrooted.pdf"
)


plot_unrooted_tree(
    standard_tree,
    standard_tree_pdf,
    "Standard amino-acid phylogeny"
)

plot_unrooted_tree(
    structural_tree,
    structural_tree_pdf,
    "3Di structural phylogeny"
)


# ============================================================
# KDE
# ============================================================

def gaussian_kde_manual(values, x_grid):
    """
    Simple Gaussian KDE using Silverman's bandwidth rule.

    Implemented directly with NumPy so scipy/seaborn are not
    required.
    """

    values = np.asarray(values, dtype=float)

    n = len(values)

    if n < 2:
        raise ValueError(
            "At least two distances are required for KDE."
        )

    std = np.std(values, ddof=1)

    # Protect against zero variance
    if std == 0:
        bandwidth = 0.5
    else:
        bandwidth = 1.06 * std * n ** (-1 / 5)

    bandwidth = max(bandwidth, 0.25)

    differences = (
        x_grid[:, None] - values[None, :]
    ) / bandwidth

    density = np.exp(
        -0.5 * differences ** 2
    ).sum(axis=1)

    density /= (
        n
        * bandwidth
        * math.sqrt(2 * math.pi)
    )

    return density


# ============================================================
# DISTANCE DISTRIBUTION PLOT
# ============================================================

standard_values = np.array([
    x["topological_distance"]
    for x in standard_distances
])

structural_values = np.array([
    x["topological_distance"]
    for x in structural_distances
])


if len(standard_values) == 0 or len(structural_values) == 0:

    raise ValueError(
        "No pairwise convergent-tip distances were generated."
    )


min_distance = min(
    standard_values.min(),
    structural_values.min()
)

max_distance = max(
    standard_values.max(),
    structural_values.max()
)


x_grid = np.linspace(
    max(0, min_distance - 2),
    max_distance + 2,
    500
)


standard_density = gaussian_kde_manual(
    standard_values,
    x_grid
)

structural_density = gaussian_kde_manual(
    structural_values,
    x_grid
)


fig, ax = plt.subplots(
    figsize=(8, 6)
)


ax.plot(
    x_grid,
    standard_density,
    linewidth=2,
    label="Standard amino-acid tree"
)

ax.plot(
    x_grid,
    structural_density,
    linewidth=2,
    label="3Di structural tree"
)


ax.set_xlabel(
    "Topological distance between convergent tips\n(number of edges)"
)

ax.set_ylabel(
    "Density"
)

ax.set_title(
    "Pairwise distances among convergent taxa"
)

ax.legend(
    frameon=False
)

plt.tight_layout()


kde_pdf = os.path.join(
    output_dir,
    "convergent_tip_distance_KDE.pdf"
)


plt.savefig(
    kde_pdf,
    format="pdf",
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# SUMMARY STATISTICS
# ============================================================

def print_summary(name, values):

    print()
    print(name)
    print("-" * len(name))

    print(f"Pairs:   {len(values)}")
    print(f"Mean:    {np.mean(values):.3f}")
    print(f"Median:  {np.median(values):.3f}")
    print(f"Min:     {np.min(values)}")
    print(f"Max:     {np.max(values)}")
    print(f"SD:      {np.std(values, ddof=1):.3f}")


print()
print("============================================================")
print("Pairwise topological-distance summary")
print("============================================================")

print_summary(
    "Standard tree",
    standard_values
)

print_summary(
    "3Di structural tree",
    structural_values
)


# ============================================================
# OUTPUT SUMMARY
# ============================================================

print()
print("============================================================")
print("Analysis completed successfully")
print("============================================================")
print()

print("Tree figures:")
print(f"  {standard_tree_pdf}")
print(f"  {structural_tree_pdf}")

print()
print("Distance tables:")
print(f"  {standard_distance_file}")
print(f"  {structural_distance_file}")

print()
print("Distance KDE:")
print(f"  {kde_pdf}")

print()