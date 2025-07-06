# build.py

import argparse
from pdbfixer import PDBFixer
from openmm.app import Modeller, PDBFile, ForceField
from openmm.unit import nanometer
import openmm as mm
import sys

def build_system(pdb_id, output_prefix="output", box_padding=1.0, ml_potential=False, ml_model_path=None):
    # Step 1: Download and fix the PDB
    print(f"Downloading and processing PDB: {pdb_id}")
    fixer = PDBFixer(pdbid=pdb_id)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    fixer.removeHeterogens(True)  # True: remove water and all heterogens

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

    # === ML POTENTIAL INTEGRATION (TorchMD-NET) ===
    if ml_potential:
        print("ML potential option selected. Swapping bonded and van der Waals terms for TorchMD-NET model.")
        try:
            from torchmdnet.openmm import TorchMDForce
        except ImportError:
            print("ERROR: TorchMD-NET is not installed. Install it with 'pip install torchmd-net'.")
            sys.exit(1)
        if not ml_model_path:
            print("ERROR: Must provide --ml_model_path when --ml_potential is enabled.")
            sys.exit(1)
        # Remove bonded and van der Waals (LJ) terms
        forces_to_remove = []
        for i, force in enumerate(system.getForces()):
            if force.__class__.__name__ in [
                "HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce", "NonbondedForce"
            ]:
                forces_to_remove.append(i)
        for idx in sorted(forces_to_remove, reverse=True):
            system.removeForce(idx)
        # Add TorchMD-NET force
        torchmd_force = TorchMDForce(ml_model_path)
        system.addForce(torchmd_force)
        print("TorchMD-NET force successfully added.")

    # Step 5: Save outputs
    with open(f"{output_prefix}_structure.pdb", "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    with open(f"{output_prefix}_system.xml", "w") as f:
        f.write(mm.XmlSerializer.serialize(system))

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
    parser.add_argument("--pdb", help="PDB ID (e.g., 7KWK)", required=True)
    parser.add_argument("--prefix", help="Output file prefix", default="output")
    parser.add_argument("--box", type=float, help="Solvent box padding in nm", default=1.0)
    parser.add_argument("--ml_potential", action="store_true", help="Use TorchMD-NET ML potential")
    parser.add_argument("--ml_model_path", type=str, help="Path to TorchMD-NET model.pt file")
    args = parser.parse_args()
    build_system(
        args.pdb, 
        args.prefix, 
        args.box,
        ml_potential=args.ml_potential,
        ml_model_path=args.ml_model_path
    )
