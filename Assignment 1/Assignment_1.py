# Nakul Sharma - 22CS10046
# Computational Biophysics
# Assigmen

import os 
import requests as req
from Bio import PDB ,SeqIO ,Align
from Bio.Seq import Seq
# from Bio.PDB.Polypeptide import th
import subprocess as sub


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
def extract_pdb_sequences(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    chains = {}

    for chain in structure.get_chains():
        chain_id = chain.get_id()
        seq = "".join(AA_MAP.get(res.get_resname().strip(), "X") for res in chain.get_residues() if res.id[0] == " ")
        chains[chain_id] = seq

    return chains
# part 2.2 : Extract the FASTA sequnces 
def read_fasta_sequences(fasta_file):
    """Read sequences from a FASTA file."""
    fasta_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_sequences[record.id] = str(record.seq)
    return fasta_sequences

#part 2.3 : compare the sequnces 
from Bio import Align
from Bio.Align import substitution_matrices

def compare_sequences(pdb_sequences, fasta_sequences):
    """Compare PDB sequences with FASTA sequences efficiently, calculating alignment, identity, and similarity."""
    
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"  # Global alignment
    matrix = substitution_matrices.load("BLOSUM62")  # Load BLOSUM62 substitution matrix

    for pdb_chain, pdb_seq in pdb_sequences.items():
        print(f"\n🔹 Chain {pdb_chain}:")

        best_match = None
        best_score = float("-inf")
        best_alignment = None
        best_identity = 0
        best_similarity = float("-inf")

        for fasta_id, fasta_seq in fasta_sequences.items():
            alignments = aligner.align(pdb_seq, fasta_seq)  # Compute alignments
            
            if not alignments:
                continue

            alignment = alignments[0]  # Take the best alignment
            score = alignment.score
            
            # Sequence Identity Calculation
            identity = sum(a == b for a, b in zip(pdb_seq, fasta_seq)) / max(len(pdb_seq), len(fasta_seq)) * 100
            
            # Sequence Similarity Calculation using BLOSUM62
            similarity = sum(matrix[a][b] for a, b in zip(pdb_seq, fasta_seq) if a in matrix and b in matrix[a])

            # Track best match
            if score > best_score:
                best_score = score
                best_match = fasta_id
                best_alignment = alignment
                best_identity = identity
                best_similarity = similarity

        # Print results for the best match
        if best_match:
            print(f"  ✅ Best match: {best_match} (Score: {best_score:.2f})")
            print(f"  🔹 Sequence Identity: {best_identity:.2f}%")
            print(f"  🔹 Sequence Similarity (BLOSUM62): {best_similarity:.2f}")
            print(best_alignment)  # Print best alignment
        else:
            print("  ❌ No match found.")



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
    pdb_sequnces = extract_pdb_sequences(PDB_F)
    fasta_seq = read_fasta_sequences(FASTA)
    compare_sequences(pdb_sequnces,fasta_seq)
    # print(pdb_sequnces)
    # print(fasta_seq)
    print("Run Successful")
