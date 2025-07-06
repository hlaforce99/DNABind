# build.py

import argparse
from pdbfixer import PDBFixer
from openmm.app import Modeller, PDBFile, ForceField
from openmm.unit import nanometer
import openmm as mm
import sys

def build_system(pdb_id, output_prefix="output", box_padding=1.0):
    # Step 1: Download and fix the PDB
    print(f"Downloading and processing PDB: {pdb_id}")
    fixer = PDBFixer(pdbid=pdb_id)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)
    
    # Step 2: Write fixed structure to temp file
    with open(f"{output_prefix}_fixed.pdb", "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    # Step 3: Load with OpenMM Modeller and solvate
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(
        ForceField('amber14/tip3p.xml'), 
        model='tip3p', 
        boxSize=None, 
        padding=box_padding*nanometer, 
        neutralize=True
    )
    
    # Step 4: Assign force fields
    # AMBER ff14SB for protein, OL15 for DNA/RNA
    forcefield = ForceField(
        'amber14-all.xml',    # ff14SB for protein, OL15 for nucleic acids
        'amber14/tip3p.xml'
    )
    print("Building system with AMBER force fields...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=mm.app.PME,
        nonbondedCutoff=1.0*nanometer,
        constraints=mm.app.HBonds
    )
    
    # Step 5: Save outputs
    with open(f"{output_prefix}_structure.pdb", "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    with open(f"{output_prefix}_system.xml", "w") as f:
        mm.XmlSerializer.serialize(system, f)
    
    # Save a default integrator for convenience
    integrator = mm.LangevinIntegrator(
        300,        # temperature in Kelvin
        1.0,        # friction coefficient in 1/ps
        0.002       # timestep in ps
    )
    with open(f"{output_prefix}_integrator.xml", "w") as f:
        f.write(mm.XmlSerializer.serialize(integrator))
    
    print(f"System built and saved as {output_prefix}_structure.pdb and {output_prefix}_system.xml")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build an OpenMM-ready system from a PDB ID.")
    parser.add_argument("--pdb", help="PDB ID (e.g., 1D66)", required=True)
    parser.add_argument("--prefix", help="Output file prefix", default="output")
    parser.add_argument("--box", type=float, help="Solvent box padding in nm", default=1.0)
    args = parser.parse_args()
    build_system(args.pdb, args.prefix, args.box)

