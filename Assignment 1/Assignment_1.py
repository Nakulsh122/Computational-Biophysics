# Nakul Sharma - 22CS10046
# Computational Biophysics
# Assigmen

import sys
import requests as req
from Bio import PDB ,SeqIO ,Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
# from Bio.PDB.Polypeptide import th
import subprocess as sub

output_file = open("report.txt","w+")
sys.stdout = output_file
AA_MOLECULAR_WEIGHTS = {
    'A': 89.094,    # Alanine
    'R': 174.203,   # Arginine
    'N': 132.119,   # Asparagine
    'D': 133.104,   # Aspartic acid
    'C': 121.154,   # Cysteine
    'E': 147.131,   # Glutamic acid
    'Q': 146.146,   # Glutamine
    'G': 75.067,    # Glycine
    'H': 155.156,   # Histidine
    'I': 131.175,   # Isoleucine
    'L': 131.175,   # Leucine
    'K': 146.189,   # Lysine
    'M': 149.208,   # Methionine
    'F': 165.192,   # Phenylalanine
    'P': 115.132,   # Proline
    'S': 105.093,   # Serine
    'T': 119.119,   # Threonine
    'W': 204.228,   # Tryptophan
    'Y': 181.191,   # Tyrosine
    'V': 117.148,   # Valine
    'X': 110.0,     # Unknown (average value)
}

AA_MAP= {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M'  # Selenomethionine is treated as Methionine
}

# downloading the required files 
def download(file_id):
    pdb_file = f"{file_id}.pdb"
    base_url_pdb = f"https://files.rcsb.org/download/{pdb_file}"
    print(base_url_pdb)
    print(f"Downloading PDB file {pdb_file}")
    response = req.get(base_url_pdb)
    if response.status_code == 200:
        with open(f"{pdb_file}","wb") as file :
           file.write(response.content)
        print("The file has been downloaded !!!")
        file.close()
    else:
        print("Error locating the file.Please Recheck your file id.")
        print(response)

def download_fasta(file_id):
    fasta_file = f"{file_id}.fasta"
    base_url_fasta = f"https://www.rcsb.org/fasta/entry/{file_id}/download"
    print(f"Downloading Fasta file {fasta_file}")
    response = req.get(base_url_fasta)
    if response.status_code == 200 :
        with open(f"{fasta_file}","wb") as fasta:
            fasta.write(response.content)
        print(f"{fasta_file} successfully downloaded")
        fasta.close()
    else:
        print("Error locating the file.Unable to download")
        print(response)
# part 2.1 : Extract protein sequence from the files .
from Bio import PDB

AA_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}


def extract_pdb_sequences(pdb_file, pdb_id):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    chains = {}
    chain_breaks = {}
    unique_chains = set()

    for chain in structure.get_chains():
        chain_id = chain.get_id()
        unique_chains.add(chain_id)

        seq = []
        residue_positions = []
        
        for res in chain.get_residues():
            if res.id[0] == " ":  # Ignore heteroatoms and water
                res_name = res.get_resname().strip()
                seq.append(AA_MAP.get(res_name, "X"))  # Convert to one-letter code
                residue_positions.append(res.id[1])   # Store residue sequence number
        
        chains[chain_id] = "".join(seq)

        # Detect chain breaks
        breaks = []
        for i in range(1, len(residue_positions)):
            if residue_positions[i] != residue_positions[i - 1] + 1:
                breaks.append(residue_positions[i - 1])  # Last residue before the break
        
        if breaks:
            chain_breaks[chain_id] = breaks  # Store breaks only if found

    num_chains = len(unique_chains)

    # Print chain break results in required format
    if chain_breaks:
        for chain_id, breaks in chain_breaks.items():
            print(f"Chain Break Detected in {pdb_id}|Chains: {chain_id}, Residue positions {breaks}")
    else:
        print(f"No Chain Break Detected in {pdb_id}")

    return chains, num_chains, chain_breaks
  # Return sequences, count, and breaks



# part 2.2 : Extract the FASTA sequnces 
def read_fasta_sequences(fasta_file):
    """Read sequences from a FASTA file."""
    fasta_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_sequences[record.id] = str(record.seq)
    return fasta_sequences

#part 2.3 : compare the sequnces 

def align_sequence_clusters(fasta_clusters, pdb_clusters):

    # Initialize aligners
    global_aligner = Align.PairwiseAligner()
    global_aligner.mode = 'global'  # Needleman-Wunsch

    local_aligner = Align.PairwiseAligner()
    local_aligner.mode = 'local'  # Smith-Waterman

    alignment_results = {}

    # Iterate through each FASTA sequence cluster
    for fasta_id, fasta_seq in fasta_clusters.items():
        fasta_seq_obj = Seq(fasta_seq)
        alignment_results[fasta_id] = {}

        # Compare with each PDB sequence cluster
        for pdb_id, pdb_seq in pdb_clusters.items():
            pdb_seq_obj = Seq(pdb_seq)

            # Perform global and local alignments
            global_alignment = global_aligner.align(fasta_seq_obj, pdb_seq_obj)[0]
            local_alignment = local_aligner.align(fasta_seq_obj, pdb_seq_obj)[0]

            # Store results
            alignment_results[fasta_id][pdb_id] = {
                "global_alignment": str(global_alignment),
                "global_score": global_alignment.score,
                "local_alignment": str(local_alignment),
                "local_score": local_alignment.score
            }

    return alignment_results

# Blosum Score Calculations 
blosum_matrices = {
    "BLOSUM62": substitution_matrices.load("BLOSUM62"),
    "BLOSUM80": substitution_matrices.load("BLOSUM80"),
    "BLOSUM45": substitution_matrices.load("BLOSUM45"),
}

def calculate_blosum_scores(fasta_dict, pdb_dict):
    
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

# main function 
if __name__ == "__main__":
    pdb_id = input("Enter the PDB id for the file :").strip().upper()
    fasta_id = input("Enter The fasta id : ").strip().upper()
    if pdb_id.endswith(".pdb"):
        pdb_id = pdb_id[:-4]
    download(pdb_id)
    download_fasta(fasta_id)
    FASTA = open(f"{fasta_id}.fasta","r")
    PDB_F = open(f"{pdb_id}.pdb","r")
    pdb_sequnces,num_chains,chain_breaks = extract_pdb_sequences(PDB_F,pdb_id)
    fasta_seq = read_fasta_sequences(FASTA)
    print(pdb_sequnces ,"\n" ,fasta_seq)
    print(f"The number of chains in the pdb sequnce {pdb_id} are :{num_chains}")
    results = align_sequence_clusters(fasta_clusters=fasta_seq,pdb_clusters=pdb_sequnces)

    for fasta_id, pdb_results in results.items():
        for pdb_id, data in pdb_results.items():
            print(f"\n Alignment between FASTA {fasta_id} and PDB {pdb_id}")
            print("\n GLOBAL ALIGNMENT (Needleman-Wunsch)")
            print(f"Score: {data['global_score']}")
            print(data["global_alignment"])
            print("-" * 50)

            print("\n LOCAL ALIGNMENT (Smith-Waterman)")
            print(f"Score: {data['local_score']}")
            print(data["local_alignment"])
            print("=" * 50)

    blosum_scores = calculate_blosum_scores(pdb_dict=pdb_sequnces,fasta_dict=fasta_seq)
    
    for matrix, results in blosum_scores.items():
        print(f"\n{matrix} Scores:")
        for key, score in results.items():
            print(f"{key}: {score}")
    print(chain_breaks)
    print("Run Successful")
    sys.stdout = sys.__stdout__
    output_file.close()
    print("report is saved ")
