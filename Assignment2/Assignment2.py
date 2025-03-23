# Assignment2 _ Nakul Sharma _ 22CS10046
# in the naccess path everywhere please put in the path to you naccess file
import os 
import subprocess
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import NeighborSearch
from Bio.SVDSuperimposer import SVDSuperimposer

def parse(inputFile):
    parser = PDBParser(QUIET=True)
    proteinStructure = parser.get_structure('protein', inputFile)
    return proteinStructure

def extract_coords(proteinStructure):
    receptorCoords = []
    ligandCoords = []
    
    for model in proteinStructure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if chain.id == 'A':
                        receptorCoords.append(atom.coord)
                    elif chain.id == 'B':
                        ligandCoords.append(atom.coord)
    
    return np.array(receptorCoords), np.array(ligandCoords)

def write_chain_pdb(proteinStructure, chainId, outputFilename):
    io = PDB.PDBIO()
    class ChainSelector(PDB.Select):
        def accept_chain(self, chain):
            return chain.get_id() == chainId
    io.set_structure(proteinStructure)
    io.save(outputFilename, ChainSelector())

def run_naccess(pdbFile, naccessPath="naccess"):
    subprocess.run(f"{naccessPath} {pdbFile}", shell=True, check=True)
    asaFile = pdbFile.replace(".pdb", ".rsa")
    totalAsa = 0.0
    if os.path.exists(asaFile):
        with open(asaFile, "r") as f:
            for line in f:
                if line.startswith("TOTAL"):
                    totalAsa = float(line.split()[2])
    return totalAsa

def compute_int_area(proteinStructure, chainA, chainB, naccessPath="naccess"):
    # Write full complex PDB
    complexPdb = "complex.pdb"
    io = PDB.PDBIO()
    io.set_structure(proteinStructure)
    io.save(complexPdb)

    complexAsa = run_naccess(complexPdb, naccessPath)

    chainAPdb = f"chain_{chainA}.pdb"
    chainBPdb = f"chain_{chainB}.pdb"
    
    write_chain_pdb(proteinStructure, chainA, chainAPdb)
    write_chain_pdb(proteinStructure, chainB, chainBPdb)

    # Compute ASA for individual chains
    asaA = run_naccess(chainAPdb, naccessPath)
    asaB = run_naccess(chainBPdb, naccessPath)

    # Compute interface area
    interfaceArea = 0.5*((asaA + asaB) - complexAsa)

    # Cleanup temporary files
    for f in [complexPdb, chainAPdb, chainBPdb, "complex.rsa", f"chain_{chainA}.rsa", f"chain_{chainB}.rsa"]:
        if os.path.exists(f):
            os.remove(f)

    return interfaceArea

def compute_solvation_energy(proteinStructure):
    solvationParams = {
        'C': 0.012, 'S': -0.214, 'O': -0.577, 'N': -0.324,
        'Fe': -8.08, 'Zn': -1.70, 'Ca': -2.44, 'Mg': -3.88,
        'P': -0.303, 'F': -0.37, 'Cl': -0.15, 'Br': -0.12, 'I': -0.10
    }
    defaultAsa = {
        'C': 15.0, 'S': 20.0, 'O': 15.0, 'N': 15.0,
        'Fe': 10.0, 'Zn': 10.0, 'Ca': 10.0, 'Mg': 10.0,
        'P': 20.0, 'F': 15.0, 'Cl': 15.0, 'Br': 15.0, 'I': 15.0,
        'X': 15.0
    }
    cutoffDist = 5.0
    
    # Get atoms separated by chain
    receptorAtoms = []
    ligandAtoms = []
    
    for model in proteinStructure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if chain.id == 'A':
                        receptorAtoms.append(atom)
                    elif chain.id == 'B':
                        ligandAtoms.append(atom)
    
    nsReceptor = NeighborSearch(receptorAtoms)
    nsLigand = NeighborSearch(ligandAtoms)
    
    interfaceLigandAtoms = []
    for atom in ligandAtoms:
        nearbyReceptorAtoms = nsReceptor.search(atom.coord, cutoffDist)
        if nearbyReceptorAtoms: 
            interfaceLigandAtoms.append(atom)
    
    interfaceReceptorAtoms = []
    for atom in receptorAtoms:
        nearbyLigandAtoms = nsLigand.search(atom.coord, cutoffDist)
        if nearbyLigandAtoms:  
            interfaceReceptorAtoms.append(atom)

    interfaceAtoms = interfaceReceptorAtoms + interfaceLigandAtoms
    solvationEnergy = 0.0
    
    for atom in interfaceAtoms:
        element = atom.element
        param = solvationParams.get(element, 0.0)
        asa = defaultAsa.get(element, defaultAsa['X'])
        energy = param * asa
        solvationEnergy += energy

    return solvationEnergy

def compute_lrmsd(targetStructure, decoyStructure):
    # Extract chains from target
    targetReceptor = targetStructure[0]['A']
    targetLigand = targetStructure[0]['B']
    
    # Extract chains from decoy
    decoyReceptor = decoyStructure[0]['A']
    decoyLigand = decoyStructure[0]['B']
     
    targetReceptorCoords = np.array([atom.coord for atom in targetReceptor.get_atoms() 
                                   if atom.name == 'CA'])
    targetLigandCoords = np.array([atom.coord for atom in targetLigand.get_atoms()])
    
    decoyReceptorCoords = np.array([atom.coord for atom in decoyReceptor.get_atoms() 
                                  if atom.name == 'CA'])
    decoyLigandCoords = np.array([atom.coord for atom in decoyLigand.get_atoms()])

    # Handle empty arrays or mismatched lengths
    if (len(targetReceptorCoords) == 0 or len(decoyReceptorCoords) == 0 or
        len(targetLigandCoords) == 0 or len(decoyLigandCoords) == 0):
        return float('inf')

    minReceptorLength = min(len(targetReceptorCoords), len(decoyReceptorCoords))
    minLigandLength = min(len(targetLigandCoords), len(decoyLigandCoords))

    targetReceptorCoords = targetReceptorCoords[:minReceptorLength]
    decoyReceptorCoords = decoyReceptorCoords[:minReceptorLength]
    targetLigandCoords = targetLigandCoords[:minLigandLength]
    decoyLigandCoords = decoyLigandCoords[:minLigandLength]

    sup = SVDSuperimposer()
    sup.set(targetReceptorCoords, decoyReceptorCoords)
    sup.run()

    rotation, transform = sup.get_rotran()

    transformedLigandDecoy = np.dot(decoyLigandCoords, rotation) + transform

    sqDiff = np.sum((targetLigandCoords - transformedLigandDecoy)**2, axis=1)
    lrmsd = np.sqrt(np.mean(sqDiff))

    return lrmsd

def compute_irmsd(targetStructure, decoyStructure):
    # Extract chains from target
    targetReceptor = targetStructure[0]['A']
    targetLigand = targetStructure[0]['B']
    
    # Extract chains from decoy
    decoyReceptor = decoyStructure[0]['A']
    decoyLigand = decoyStructure[0]['B']
    
    targetReceptorAtoms = list(targetReceptor.get_atoms())
    targetLigandAtoms = list(targetLigand.get_atoms())
    decoyReceptorAtoms = list(decoyReceptor.get_atoms())
    decoyLigandAtoms = list(decoyLigand.get_atoms())
    cutoffDist = 5.0

    nsTargetReceptor = NeighborSearch(targetReceptorAtoms)
    targetInterfaceLigandResidues = set()
    for atom in targetLigandAtoms:
        closeAtoms = nsTargetReceptor.search(atom.coord, cutoffDist)
        if closeAtoms:
            targetInterfaceLigandResidues.add(atom.get_parent())
    
    nsTargetLigand = NeighborSearch(targetLigandAtoms)
    targetInterfaceReceptorResidues = set()
    for atom in targetReceptorAtoms:
        closeAtoms = nsTargetLigand.search(atom.coord, cutoffDist)
        if closeAtoms:
            targetInterfaceReceptorResidues.add(atom.get_parent())
    
    nsDecoyReceptor = NeighborSearch(decoyReceptorAtoms)
    decoyInterfaceLigandResidues = set()
    for atom in decoyLigandAtoms:
        closeAtoms = nsDecoyReceptor.search(atom.coord, cutoffDist)
        if closeAtoms:
            decoyInterfaceLigandResidues.add(atom.get_parent())
    
    nsDecoyLigand = NeighborSearch(decoyLigandAtoms)
    decoyInterfaceReceptorResidues = set()
    for atom in decoyReceptorAtoms:
        closeAtoms = nsDecoyLigand.search(atom.coord, cutoffDist)
        if closeAtoms:
            decoyInterfaceReceptorResidues.add(atom.get_parent())
    
    targetInterfaceCas = []
    for res in targetInterfaceReceptorResidues.union(targetInterfaceLigandResidues):
        for atom in res:
            if atom.name == 'CA':
                targetInterfaceCas.append(atom.coord)
    
    decoyInterfaceCas = []
    for res in decoyInterfaceReceptorResidues.union(decoyInterfaceLigandResidues):
        for atom in res:
            if atom.name == 'CA':
                decoyInterfaceCas.append(atom.coord)
    
    if not targetInterfaceCas or not decoyInterfaceCas:
        return float('inf') 
    
    minLength = min(len(targetInterfaceCas), len(decoyInterfaceCas))
    targetInterfaceCas = np.array(targetInterfaceCas[:minLength])
    decoyInterfaceCas = np.array(decoyInterfaceCas[:minLength])
    
    sup = SVDSuperimposer()
    sup.set(targetInterfaceCas, decoyInterfaceCas)
    sup.run()
    
    return sup.get_rms()

def compute_fnat(targetStructure, decoyStructure, cutoffDist=5.0):
    # Extract chains from target
    targetReceptor = targetStructure[0]['A']
    targetLigand = targetStructure[0]['B']
    
    # Extract chains from decoy
    decoyReceptor = decoyStructure[0]['A']
    decoyLigand = decoyStructure[0]['B']
    
    targetReceptorAtoms = list(targetReceptor.get_atoms())
    targetLigandAtoms = list(targetLigand.get_atoms())
    decoyReceptorAtoms = list(decoyReceptor.get_atoms())
    decoyLigandAtoms = list(decoyLigand.get_atoms())
    
    nsTargetReceptor = NeighborSearch(targetReceptorAtoms)
    nativeContacts = set()
    
    for ligandAtom in targetLigandAtoms:
        closeReceptorAtoms = nsTargetReceptor.search(ligandAtom.coord, cutoffDist)
        for receptorAtom in closeReceptorAtoms:
            contact = (
                receptorAtom.get_parent().get_id(), 
                ligandAtom.get_parent().get_id()
            )
            nativeContacts.add(contact)
    
    nsDecoyReceptor = NeighborSearch(decoyReceptorAtoms)
    decoyContacts = set()
    
    for ligandAtom in decoyLigandAtoms:
        closeReceptorAtoms = nsDecoyReceptor.search(ligandAtom.coord, cutoffDist)
        for receptorAtom in closeReceptorAtoms:
            contact = (
                receptorAtom.get_parent().get_id(), 
                ligandAtom.get_parent().get_id()
            )
            decoyContacts.add(contact)
    
    if len(nativeContacts) == 0:
        return 0.0
    
    commonContacts = nativeContacts.intersection(decoyContacts)
    
    fnat = len(commonContacts) / len(nativeContacts)
    
    return fnat

def evaluate(targetFile, folderPath, reportFile):
    # Parse target structure
    targetStructure = parse(targetFile)
    
    with open(reportFile, 'w') as f:
        f.write("Filename\t | Interface Area (A^2)\t | Solvation Energy (kcal/mol)\t | LRMSD (A)\t | IRMSD (A)\t | Fnat\n")
        f.write("-"*80 + "\n")

        for decoyFilename in os.listdir(folderPath):
            decoyPath = os.path.join(folderPath, decoyFilename)
            
            # Parse decoy structure
            decoyStructure = parse(decoyPath)

            interfaceArea = compute_int_area(decoyStructure, "A", "B")
            solvationEnergy = compute_solvation_energy(decoyStructure)
            lrmsdVal = compute_lrmsd(targetStructure, decoyStructure)
            irmsdVal = compute_irmsd(targetStructure, decoyStructure)
            fnatVal = compute_fnat(targetStructure, decoyStructure)

            f.write(f"{decoyFilename}\t | {interfaceArea:.2f}\t | {solvationEnergy:.2f}\t | {lrmsdVal:.2f}\t | {irmsdVal:.2f}\t | {fnatVal:.2f}\n")
            f.write("-"*80 + "\n")

if __name__ == "__main__":
    targetPdb = "target.pdb"
    decoyFolder = "Decoys"
    reportFile = "Score.txt"

    evaluate(targetPdb, decoyFolder, reportFile)
    print("Run Successful, file saved")