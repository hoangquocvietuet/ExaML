import os
import sys
import subprocess
from Bio import SeqIO

def generate_control_file(sites, taxa, partitions, output_name):
    """
    Generates an INDELible control file to ensure maximum uniqueness.
    """
    control_filename = f"control.txt"
    phylip_filename = f"{output_name}_TRUE.phy"
    partition_filename = f"{output_name}_partitions.txt"

    with open(control_filename, "w") as f:
        f.write("/////////////////////////////////////////////////////////////////////////////////////\n")
        f.write("// INDELible V1.03 control file - Auto-generated for Maximum Uniqueness\n")
        f.write(f"// Simulating {taxa} taxa, {sites} sites, {partitions} partitions\n")
        f.write("/////////////////////////////////////////////////////////////////////////////////////\n\n")

        f.write("[TYPE] NUCLEOTIDE 1\n\n")

        f.write("[SETTINGS]\n")
        f.write("    [output] PHYLIP\n")
        f.write("    [randomseed] 4321\n\n")

        # Define highly diverse models for each partition
        for i in range(1, partitions + 1):
            f.write(f"[MODEL] Model{i}\n")
            f.write("    [submodel] GTR 10.0 10.0 10.0 10.0 10.0 10.0\n")  # Extreme mutation bias
            f.write("    [statefreq] 0.25 0.25 0.25 0.25\n")
            f.write("    [rates] 0.0 10.0 4\n")  # Extreme rate heterogeneity
            f.write("    [insertrate] 0.0\n")
            f.write("    [deleterate] 0.0\n\n")

        # Define an extremely diverse tree
        f.write("[TREE] SimTree\n")
        f.write(f"    [unrooted] {taxa} 10.0 2.0 1.0 2.0\n")  # High birth-death process
        f.write("    [treelength] 7.0\n\n")  # Long tree = many mutations

        # Define partitions
        f.write("[PARTITIONS] AutoPartition\n")
        partition_sizes = sites // partitions
        for i in range(1, partitions + 1):
            f.write(f"    [SimTree Model{i} {partition_sizes}]\n")
        f.write("\n")

        # Run evolution
        f.write("[EVOLVE]\n")
        f.write(f"    AutoPartition 1 {output_name}\n")

    print(f"Generated control file: {control_filename}")
    return control_filename, phylip_filename, partition_filename, partition_sizes

def run_indelible(output_name):
    """
    Runs INDELible using the generated control file.
    """
    print(f"Running INDELible for {output_name}...")
    log_file = f"{output_name}.log"
    with open(log_file, "w") as log:
        subprocess.run(["./indelible"], stdout=log, stderr=log)
    print(f"INDELible simulation completed for {output_name}.")

def read_fasta(fasta_file):
    """
    Reads the generated FASTA file and calculates:
    - Number of unique site patterns
    - Number of unique sequences
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))  # Convert Seq object to string

    if len(sequences) == 0:
        print("No sequences found in the output file.")
        return 0, 0

    # Check sequence length consistency
    seq_length = len(sequences[0])
    for seq in sequences:
        if len(seq) != seq_length:
            print(f"Sequences have different lengths: {len(seq)} vs {seq_length}")
            return 0, 0

    # Extract columns
    columns = [[] for _ in range(seq_length)]
    for seq in sequences:
        for i in range(seq_length):
            columns[i].append(seq[i])

    # Count unique patterns (site columns)
    pattern = set(tuple(col) for col in columns)

    # Count distinct sequences
    distinct_seqs = set(sequences)

    return len(pattern), len(distinct_seqs)

def read_phylip(phylip_file):
    """
    Reads a PHYLIP file and calculates:
    - Number of unique site patterns
    - Number of unique sequences
    """
    sequences = []
    for record in SeqIO.parse(phylip_file, "phylip"):
        sequences.append(str(record.seq))  # Convert Seq object to string

    if len(sequences) == 0:
        print("No sequences found in the PHYLIP file.")
        return 0, 0

    # Check sequence length consistency
    seq_length = len(sequences[0])
    for seq in sequences:
        if len(seq) != seq_length:
            print(f"Sequences have different lengths: {len(seq)} vs {seq_length}")
            return 0, 0

    # Extract columns
    columns = [[] for _ in range(seq_length)]
    for seq in sequences:
        for i in range(seq_length):
            columns[i].append(seq[i])

    # Count unique patterns (site columns)
    pattern = set(tuple(col) for col in columns)

    # Count distinct sequences
    distinct_seqs = set(sequences)

    return len(pattern), len(distinct_seqs)

def write_partition_file(sites, partitions, partition_filename):
    """
    Writes the partition file in the format:
    DNA, part1 = 1-100
    DNA, part2 = 101-384
    """
    partition_sizes = sites // partitions
    with open(partition_filename, "w") as f:
        for i in range(1, partitions + 1):
            start = (i - 1) * partition_sizes + 1
            end = start + partition_sizes - 1
            f.write(f"DNA, part{i} = {start}-{end}\n")
    
    print(f"Partition file saved as: {partition_filename}")

def move_output_to_folder(output_name):
    """
    Moves generated output files to a folder named after the output_name.
    """
    folder_name = f"./{output_name}_results"
    os.makedirs(folder_name, exist_ok=True)

    files_to_move = [
        f"{output_name}_TRUE.phy",
        f"{output_name}_partitions.txt",
        f"control.txt",
        f"{output_name}.log",
        f"newick_trees.txt"
    ]

    for file in files_to_move:
        if os.path.exists(file):
            os.rename(file, os.path.join(folder_name, file))

    print(f"All output files moved to {folder_name}")
    
def extract_newick_trees(trees_file="trees.txt", output_file="newick_trees.txt"):
    """
    Extracts Newick tree strings from the trees.txt file and saves them into newick_trees.txt.
    """
    newick_trees = []

    with open(trees_file, "r") as infile:
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) > 7 and "(" in parts[-1] and ");" in parts[-1]:  # Newick trees start with '(' and end with ');'
                newick_trees.append(parts[-1])  # Extract Newick tree

    # Save extracted Newick trees to a file
    with open(output_file, "w") as outfile:
        for tree in newick_trees:
            outfile.write(tree + "\n")

    print(f"✅ Extracted {len(newick_trees)} Newick trees and saved to {output_file}")

def batch_run(configurations):
    """
    Runs multiple simulations based on an array of configurations.
    Each configuration is a tuple: (sites, taxa, patterns, partitions, output_name)
    """
    for config in configurations:
        sites, taxa, partitions, output_name = config
        print(f"\n--- Running Simulation: {output_name} ---")

        # Generate control file
        control_file, phylip_output, partition_file, partition_sizes = generate_control_file(
            sites, taxa, partitions, output_name
        )

        # Run INDELible
        run_indelible(output_name)

        # Analyze the generated PHYLIP file
        num_patterns, num_unique_sequences = read_phylip(phylip_output)

        # Write partition file
        write_partition_file(sites, partitions, partition_file)
        
        extract_newick_trees()

        print(f"\nAnalysis Results for {output_name}:")
        print(f"  - Unique site patterns: {num_patterns}")
        print(f"  - Unique sequences: {num_unique_sequences}")

        # Move outputs to a separate folder
        move_output_to_folder(output_name)

if __name__ == "__main__":
    # Define multiple configurations as (sites, taxa, patterns, partitions, output_name)
    configurations = [
        (50, 5, 1, "test1"),
        # (300000, 200, 3, "test2"),
        # (300000, 200, 50, "test3"),
        # (4000000, 20, 3, "test4"),
        # (4000000, 20, 500, "test4"),
    ]

    batch_run(configurations)
