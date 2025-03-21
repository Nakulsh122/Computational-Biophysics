import os 
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import NeighborSearch
from Bio.SVDSuperimposer import SVDSuperimposer

def parse(file):
    parser = PDBParser(QUIET=True)
    structure = parser('protien',file)

    receptor = []
    ligand = []

    for modal in structure:
        for chain in modal:
            for residue in chain:
                for atom in residue:
                    if chain.id == 'A':
                        receptor.append(atom.coord)
                    elif chain.id == 'B':
                        ligand.append(atom.coord)
    return np.array(receptor),np.array(ligand)
def compute_int_area(receptor,ligand):
    cut_off = 5.0
    ns_receptors = NeighborSearch(receptor)
    interface_atoms = set()

    for ligand_atom in ligand : 
        nearby_atoms = ns_receptors.search(ligand_atom.coord,cut_off)
        if nearby_atoms :
            interface_atoms.add(ligand_atom.serial_number)
            for receptor_atom in nearby_atoms:
                interface_atoms.add(receptor_atom.serial_number)

    esitmated_area = len(interface_atoms)*10
    return esitmated_area
def compute_solvation_energy(receptor,ligand):
    solvation_params = {
        'C': 0.012, 'S': -0.214, 'O': -0.577, 'N': -0.324,
        'Fe': -8.08, 'Zn': -1.70, 'Ca': -2.44, 'Mg': -3.88,
        'P': -0.303, 'F': -0.37, 'Cl': -0.15, 'Br': -0.12, 'I': -0.10
    }
    default_asa = {
        'C': 15.0, 'S': 20.0, 'O': 15.0, 'N': 15.0,
        'Fe': 10.0, 'Zn': 10.0, 'Ca': 10.0, 'Mg': 10.0,
        'P': 20.0, 'F': 15.0, 'Cl': 15.0, 'Br': 15.0, 'I': 15.0,
        'X': 15.0
    }
    cutoff = 0.5
    receptor_atoms = list(receptor)
    ligand_atoms = list(ligand)
    ns_receptor = NeighborSearch(receptor_atoms)
    ns_ligand = NeighborSearch(ligand_atoms)
    interface_ligand_atoms = []
    for atom in ligand_atoms:
        nearby_receptor_atoms = ns_receptor.search(atom.coord, cutoff)
        if nearby_receptor_atoms: 
            interface_ligand_atoms.append(atom)
    interface_receptor_atoms = []
    for atom in receptor_atoms:
        nearby_ligand_atoms = ns_ligand.search(atom.coord, cutoff)
        if nearby_ligand_atoms:  
            interface_receptor_atoms.append(atom)

    interface_atoms = interface_receptor_atoms + interface_ligand_atoms
    solvation_energy  = 0.0
    for atom in interface_atoms :
        element=atom.element
        param = solvation_params.get(element,0.0)
        asa = default_asa.get(element,default_asa['X'])
        energy = param*asa
        solvation_energy += energy

    return solvation_energy
    pass
def compute_lrmsd(target_receptor,target_ligand,decoy_receptor,decoy_ligand):
     
    target_receptor_coords = np.array([atom.coord for atom in target_receptor.get_atoms() 
                                       if atom.name == 'CA']) 
    target_ligand_coords = np.array([atom.coord for atom in target_ligand.get_atoms()])
    
    decoy_receptor_coords = np.array([atom.coord for atom in decoy_receptor.get_atoms() 
                                      if atom.name == 'CA'])
    decoy_ligand_coords = np.array([atom.coord for atom in decoy_ligand.get_atoms()]) 

    min_receptor_length = min(len(target_receptor_coords), len(decoy_receptor_coords))
    min_ligand_length = min(len(target_ligand_coords), len(decoy_ligand_coords))

    target_receptor_coords = target_receptor_coords[:min_receptor_length]
    decoy_receptor_coords = decoy_receptor_coords[:min_receptor_length]
    target_ligand_coords = target_ligand_coords[:min_ligand_length]
    decoy_ligand_coords = decoy_ligand_coords[:min_ligand_length]

    sup = SVDSuperimposer()
    sup.set(target_receptor_coords, decoy_receptor_coords)
    sup.run()

    rotation,transform = sup.get_rotran()

    transfromed_ligand_decoy = np.dot(decoy_ligand_coords,rotation) + transform

    sq_diff = np.sum((target_ligand_coords - transfromed_ligand_decoy)**2 , axis=1)
    lrmsd = np.sqrt(np.mean(sq_diff))

    return lrmsd
def compute_irmsd(target_receptor,target_ligand,decoy_receptor,decoy_ligand):
    target_receptor_atoms = list(target_receptor.get_atoms())
    target_ligand_atoms = list(target_ligand.get_atoms())
    decoy_receptor_atoms = list(decoy_receptor.get_atoms())
    decoy_ligand_atoms = list(decoy_ligand.get_atoms())
    cutoff=0.5

    ns_target_receptor = NeighborSearch(target_receptor_atoms)
    target_interface_ligand_residues = set()
    for atom in target_ligand_atoms:
        close_atoms = ns_target_receptor.search(atom.coord, cutoff)
        if close_atoms:
            target_interface_ligand_residues.add(atom.get_parent())
    
    ns_target_ligand = NeighborSearch(target_ligand_atoms)
    target_interface_receptor_residues = set()
    for atom in target_receptor_atoms:
        close_atoms = ns_target_ligand.search(atom.coord, cutoff)
        if close_atoms:
            target_interface_receptor_residues.add(atom.get_parent())
    
    ns_decoy_receptor = NeighborSearch(decoy_receptor_atoms)
    decoy_interface_ligand_residues = set()
    for atom in decoy_ligand_atoms:
        close_atoms = ns_decoy_receptor.search(atom.coord, cutoff)
        if close_atoms:
            decoy_interface_ligand_residues.add(atom.get_parent())
    
    ns_decoy_ligand = NeighborSearch(decoy_ligand_atoms)
    decoy_interface_receptor_residues = set()
    for atom in decoy_receptor_atoms:
        close_atoms = ns_decoy_ligand.search(atom.coord, cutoff)
        if close_atoms:
            decoy_interface_receptor_residues.add(atom.get_parent())
    

    target_interface_cas = []
    for res in target_interface_receptor_residues.union(target_interface_ligand_residues):
        for atom in res:
            if atom.name == 'CA':
                target_interface_cas.append(atom.coord)
    
    decoy_interface_cas = []
    for res in decoy_interface_receptor_residues.union(decoy_interface_ligand_residues):
        for atom in res:
            if atom.name == 'CA':
                decoy_interface_cas.append(atom.coord)
    

    if not target_interface_cas or not decoy_interface_cas:
        return float('inf') 
    

    min_length = min(len(target_interface_cas), len(decoy_interface_cas))
    target_interface_cas = np.array(target_interface_cas[:min_length])
    decoy_interface_cas = np.array(decoy_interface_cas[:min_length])
    

    sup = SVDSuperimposer()
    sup.set(target_interface_cas, decoy_interface_cas)
    sup.run()
    

    return sup.get_rms()

def compute_fnat(target_receptor, target_ligand, decoy_receptor, decoy_ligand, cutoff=5.0):
    from Bio.PDB import NeighborSearch
    

    target_receptor_atoms = list(target_receptor.get_atoms())
    target_ligand_atoms = list(target_ligand.get_atoms())
    

    decoy_receptor_atoms = list(decoy_receptor.get_atoms())
    decoy_ligand_atoms = list(decoy_ligand.get_atoms())
    

    ns_target_receptor = NeighborSearch(target_receptor_atoms)
    native_contacts = set()
    
    for ligand_atom in target_ligand_atoms:
        close_receptor_atoms = ns_target_receptor.search(ligand_atom.coord, cutoff)
        for receptor_atom in close_receptor_atoms:

            contact = (
                receptor_atom.get_parent().get_id(), 
                ligand_atom.get_parent().get_id()
            )
            native_contacts.add(contact)
    

    ns_decoy_receptor = NeighborSearch(decoy_receptor_atoms)
    decoy_contacts = set()
    
    for ligand_atom in decoy_ligand_atoms:
        close_receptor_atoms = ns_decoy_receptor.search(ligand_atom.coord, cutoff)
        for receptor_atom in close_receptor_atoms:

            contact = (
                receptor_atom.get_parent().get_id(), 
                ligand_atom.get_parent().get_id()
            )
            decoy_contacts.add(contact)
    

    if len(native_contacts) == 0:
        return 0.0
    

    common_contacts = native_contacts.intersection(decoy_contacts)
    

    fnat = len(common_contacts) / len(native_contacts)
    
    return fnat

def evaluate(target,folder,file):

    target_receptor,target_ligand = parse(target)
    with open(file,'w') as f:
        f.write("Filename | Interface Area (Å²) | Solvation Energy (kcal/mol) | LRMSD (Å) | IRMSD (Å) | Fnat\n")
        f.write("-"*80 + "\n")

        for decoy in os.listdir(folder):
            path = os.path.join(folder,decoy)
            decoy_receptor,decoy_ligand = parse(path)

            int_area  = compute_int_area(decoy_receptor,decoy_ligand)
            solvation_energy = compute_solvation_energy(decoy_receptor,decoy_ligand)
            lrmsd = compute_lrmsd(target_receptor, target_ligand, decoy_receptor, decoy_ligand)
            irmsd = compute_irmsd(target_receptor, target_ligand, decoy_receptor, decoy_ligand)
            fnat = compute_fnat(target_receptor, target_ligand, decoy_receptor, decoy_ligand)

            f.write(f"{decoy} | {int_area:.2f} | {solvation_energy:.2f} | {lrmsd:.2f} | {irmsd:.2f} | {fnat:.2f}\n")
            f.write("-"*80 + "\n")
    pass

if __name__ == "__main__" :

    target_pdb = input("Enter the name of the targe pdb file:")
    Decoy_folder = input("Enter the name of the edcoyy folder:")
    report_file = "Report.txt"

    evaluate(target_pdb,Decoy_folder,report_file)
    print("Run Successful,file saved")
