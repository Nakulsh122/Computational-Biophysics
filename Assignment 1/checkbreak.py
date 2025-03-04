def detect_chain_breaks(fasta_clusters, pdb_clusters):
    """
    Detects chain breaks by comparing PDB sequences with FASTA sequences.
    """
    chain_breaks = {}

    for chain_id, fasta_seq in fasta_clusters.items():
        # Find corresponding PDB sequence
        pdb_seq = None
        for key, sequence in pdb_clusters.items():
            if f"|Chains" in key and chain_id in key:
                pdb_seq = sequence
                break

        if not pdb_seq:
            print(f"⚠️ No PDB sequence found for Chain {chain_id}, skipping...")
            continue

        # Detect chain breaks by checking missing segments
        missing_positions = []
        fasta_index = 0

        for i, res in enumerate(pdb_seq):
            if fasta_index >= len(fasta_seq):
                break  # Prevent out-of-bounds error

            if res != fasta_seq[fasta_index]:  # If mismatch, look for missing residues
                missing_positions.append(fasta_index + 1)  # 1-based indexing
                continue

            fasta_index += 1  # Move to the next FASTA residue

        if missing_positions:
            chain_breaks[chain_id] = missing_positions

    return chain_breaks

fasta_clusters = {
    'A': 'VLSPADKTNVKAAW',
    'B': 'VHLTPEEKSAVTAL'
}

pdb_clusters = {
    '1XYZ_A|Chains': 'VLSPADKTNVKAAW',
    '1XYZ_B|Chains': 'VHLTPEEKSAVTAL'
}

# Expected Output:
# "No chain breaks detected."

# Detect chain breaks
chain_breaks = detect_chain_breaks(fasta_clusters, pdb_clusters)

# Save results to report.txt
with open("report.txt", "w") as report:
    report.write("\n=== Chain Break Detection ===\n")
    if chain_breaks:
        for chain, positions in chain_breaks.items():
            report.write(f"Chain {chain} has breaks at positions: {positions}\n")
    else:
        report.write("No chain breaks detected.\n")
report.close()
print("✅ Chain break detection complete! Results saved in report.txt")
