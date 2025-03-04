# Nakul Sharma - 22CS10046
# Computational Biophysics
# Assigmen

import os 
import requests as req


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
    else:
        print("Error locating the file.Unable to download")
        print(response)

if __name__ == "__main__":
    pdb_id = input("Enter the PDB id for the file :").strip().upper()
    if pdb_id.endswith(".pdb"):
        pdb_id = pdb_id[:-4]
    download(pdb_id)
    download_fasta(pdb_id)
