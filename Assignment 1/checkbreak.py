def detect_chain_breaks(pdb_clusters, fasta_clusters):
    chain_breaks = {}

    for chain_id, fasta_seq in fasta_clusters.items():
        # Try to find a matching chain in PDB data
        matching_pdb_chain = None
        for pdb_chain in pdb_clusters:
            if chain_id in pdb_chain or pdb_chain in chain_id:
                matching_pdb_chain = pdb_chain
                break

        if not matching_pdb_chain:
            print(f"Warning: No matching PDB chain found for {chain_id}")
            continue

        pdb_seq = pdb_clusters[matching_pdb_chain]
        break_positions = []
        pdb_index = 0  # Track position in PDB sequence

        for fasta_index, fasta_residue in enumerate(fasta_seq):
            # If the PDB sequence has more residues to compare
            if pdb_index < len(pdb_seq):
                if pdb_seq[pdb_index] == fasta_residue:
                    pdb_index += 1  # Move to next PDB residue
                else:
                    break_positions.append(fasta_index + 1)  # 1-based index
            else:
                break_positions.append(fasta_index + 1)

        if break_positions:
            chain_breaks[chain_id] = break_positions    

    return chain_breaks

# Run the function

 # FASTA Clusters (Corrected)
fasta_clusters = {
    '1AON_1|Chains': 'AAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTTTATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLLIIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEYGNMIDMGILDPTKVTRSALQYAASVAGLMITTECMVTDLPKNDAADLGAAGGMGGMGGMGGMM',
    
    '1AON_2|Chains': 'MNIRPLHDRVIVKRKEVETKSAGGIVLTGSAAAKSTRGEVLAVGNGRILENGEVKPLDVKVGDIVIFNDGYGVKSEKIDNEEVLIMSESDILAIVEA'
}

# PDB Clusters (Corrected)
pdb_clusters = {
    'A': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
    'B': 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
    'C': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
    'D': 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
}

# Run the function
detect_chain_breaks(pdb_clusters, fasta_clusters)
