from Bio.Align import substitution_matrices

# Load BLOSUM matrices
blosum_matrices = {
    "BLOSUM62": substitution_matrices.load("BLOSUM62"),
    "BLOSUM80": substitution_matrices.load("BLOSUM80"),
    "BLOSUM45": substitution_matrices.load("BLOSUM45"),
}

def calculate_blosum_scores(fasta_dict, pdb_dict):
    """
    Compute BLOSUM45, BLOSUM62, and BLOSUM80 scores for each pair between FASTA and PDB clusters.
    
    Args:
        fasta_dict (dict): Dictionary of FASTA sequences with labels.
        pdb_dict (dict): Dictionary of PDB sequences with labels.

    Returns:
        dict: Nested dictionary with BLOSUM scores for each matrix and sequence pair.
    """
    scores = {matrix_name: {} for matrix_name in blosum_matrices}

    for fasta_label, fasta_seq in fasta_dict.items():
        for pdb_label, pdb_seq in pdb_dict.items():
            pair_key = (fasta_label, pdb_label)  # Maintain original labels
            min_length = min(len(fasta_seq), len(pdb_seq))

            for matrix_name, blosum_matrix in blosum_matrices.items():
                score = sum(
                    blosum_matrix.get((fasta_seq[i], pdb_seq[i]), -4)  # Default mismatch penalty
                    for i in range(min_length)
                )
                scores[matrix_name][pair_key] = score
    
    return scores
fasta_sequences = {
    'A': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
    'B': 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
    'C': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
    'D': 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
}

pdb_sequences = {
    '1AON_1|Chains': 'AAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTTTATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLLIIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEYGNMIDMGILDPTKVTRSALQYAASVAGLMITTECMVTDLPKNDAADLGAAGGMGGMGGMGGMM',
    
    '1AON_2|Chains': 'MNIRPLHDRVIVKRKEVETKSAGGIVLTGSAAAKSTRGEVLAVGNGRILENGEVKPLDVKVGDIVIFNDGYGVKSEKIDNEEVLIMSESDILAIVEA'
}



blosum_scores = calculate_blosum_scores(fasta_sequences, pdb_sequences)

for matrix, results in blosum_scores.items():
    print(f"\n{matrix} Scores:")
    for key, score in results.items():
        print(f"{key}: {score}")
