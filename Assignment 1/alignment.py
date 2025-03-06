from Bio import Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices

def align_sequences(seq1, seq2):
    """
    Aligns two sequences using both global (Needleman-Wunsch) and local (Smith-Waterman) alignment.
    Uses BLOSUM62 substitution matrix.

    Parameters:
    seq1 (str): First sequence.
    seq2 (str): Second sequence.

    Returns:
    dict: Contains global and local alignment results, scores, and similarity percentages.
    """

    # Initialize aligners
    global_aligner = Align.PairwiseAligner()
    local_aligner = Align.PairwiseAligner()

    # Set modes
    global_aligner.mode = 'global'  # Needleman-Wunsch
    local_aligner.mode = 'local'  # Smith-Waterman

    # Load BLOSUM62 matrix
    blosum62 = substitution_matrices.load("BLOSUM62")
    global_aligner.substitution_matrix = blosum62
    local_aligner.substitution_matrix = blosum62

    # Perform alignments
    global_alignment = global_aligner.align(seq1, seq2)[0]
    local_alignment = local_aligner.align(seq1, seq2)[0]

    # Compute similarity
    global_similarity = calculate_similarity(global_alignment)
    local_similarity = calculate_similarity(local_alignment)

    return {
        "global_alignment": str(global_alignment),
        "global_score": global_alignment.score,
        "global_similarity": global_similarity,
        "local_alignment": str(local_alignment),
        "local_score": local_alignment.score,
        "local_similarity": local_similarity
    }

def calculate_similarity(alignment):
    """
    Computes the similarity percentage between two aligned sequences.

    Parameters:
    alignment (Bio.Align.Alignment): The alignment result.

    Returns:
    float: Similarity percentage.
    """
    aligned_seq1 = alignment.aligned[0]
    aligned_seq2 = alignment.aligned[1]

    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    total_length = alignment.shape[1]  # Length of alignment

    return (matches / total_length) * 100 if total_length > 0 else 0

def align_sequence_clusters(fasta_clusters, pdb_clusters):
    """
    Aligns multiple clusters of FASTA sequences to PDB sequences.

    Parameters:
    fasta_clusters (dict): Dictionary of FASTA sequences {cluster_id: sequence}
    pdb_clusters (dict): Dictionary of PDB sequences {cluster_id: sequence}

    Returns:
    dict: A dictionary with alignment results.
    """

    alignment_results = {}

    for fasta_id, fasta_seq in fasta_clusters.items():
        fasta_seq_obj = Seq(fasta_seq)
        alignment_results[fasta_id] = {}

        for pdb_id, pdb_seq in pdb_clusters.items():
            pdb_seq_obj = Seq(pdb_seq)

            # Perform alignment for this pair
            alignment_data = align_sequences(fasta_seq_obj, pdb_seq_obj)
            alignment_results[fasta_id][pdb_id] = alignment_data

    return alignment_results

def format_alignment_results(results):
    """
    Formats and prints the alignment results in a readable way.

    Parameters:
    results (dict): The output from align_sequence_clusters.
    """

    for fasta_id, pdb_results in results.items():
        for pdb_id, data in pdb_results.items():
            print(f"\n🔹 Alignment between FASTA {fasta_id} and PDB {pdb_id} 🔹")

            print("\n✅ GLOBAL ALIGNMENT (Needleman-Wunsch)")
            print(f"Score: {data['global_score']}")
            print(f"Similarity: {data['global_similarity']:.2f}%")
            print(data["global_alignment"])
            print("-" * 50)

            print("\n✅ LOCAL ALIGNMENT (Smith-Waterman)")
            print(f"Score: {data['local_score']}")
            print(f"Similarity: {data['local_similarity']:.2f}%")
            print(data["local_alignment"])
            print("=" * 50)

# Example Data
fasta_clusters = {
    'A': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
    'B': 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
}

pdb_clusters = {
    '1AON_1': 'AAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTTTATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLLIIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEYGNMIDMGILDPTKVTRSALQYAASVAGLMITTECMVTDLPKNDAADLGAAGGMGGMGGMGGMM',
    '1AON_2': 'MNIRPLHDRVIVKRKEVETKSAGGIVLTGSAAAKSTRGEVLAVGNGRILENGEVKPLDVKVGDIVIFNDGYGVKSEKIDNEEVLIMSESDILAIVEA'
}

# Run the alignment function
results = align_sequence_clusters(fasta_clusters, pdb_clusters)

# Print formatted results
format_alignment_results(results)
