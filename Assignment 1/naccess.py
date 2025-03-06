import subprocess

def run_naccess(pdb_filename):
    
    naccess_path = "/mnt/d/NACCESS/naccess"  # Full path to NACCESS executable
    subprocess.run(f"{naccess_path} {pdb_filename}", shell=True)
