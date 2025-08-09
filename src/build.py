#!/usr/bin/env python3
"""
build.py

System preparation module for DNABind pipeline:
- Downloads PDB structures or uses local files
- Extracts and parameterizes ligands using AmberTools (GAFF2)
- Prepares DNA structure with missing atoms/hydrogens
- Builds and solvates the complex using tleap
- Generates OpenMM system XML files
- Optional ML potential integration (TorchMD-NET)
"""

import os
import subprocess
import warnings
from pdbfixer import PDBFixer
import openmm.app as app
import openmm as mm
import openmm.unit as unit
from utils import parse_chain_resid

# Suppress noisy warnings from OpenMM and PDBFixer
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Periodic table for electron counting
PERIODIC_TABLE = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "P": 15,
    "S": 16,
    "F": 9,
    "CL": 17,
    "BR": 35,
    "I": 53,
}


def run_command(cmd, silent=True):
    """Run a shell command and handle errors"""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
            text=True,
        )
        if not silent:
            print(result.stdout)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {cmd}")
        print(e.stdout)
        raise RuntimeError(f"Command failed with exit code {e.returncode}")


def download_pdb(pdb_id, output_file=None):
    """Download PDB file from RCSB"""
    if output_file is None:
        output_file = f"{pdb_id}.pdb"

    if os.path.exists(output_file):
        print(f"Using existing PDB file: {output_file}")
        return output_file

    print(f"Downloading PDB {pdb_id}...")
    from urllib.request import urlopen

    try:
        with urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb") as response:
            pdb_content = response.read().decode()

        with open(output_file, "w") as f:
            f.write(pdb_content)

        print(f"Downloaded PDB to {output_file}")
        return output_file
    except Exception as e:
        print(f"Error downloading PDB {pdb_id}: {e}")
        raise


def extract_ligand(pdb_file, ligand_resname=None, output_prefix="output"):
    """Extract ligand from PDB file and prepare it for parameterization"""
    standard_residues = {
        "DA",
        "DT",
        "DG",
        "DC",
        "DU",
        "A",
        "T",
        "G",
        "C",
        "U",
        "HOH",
        "WAT",
    }

    # If no ligand name provided, find the largest non-standard residue
    if ligand_resname is None:
        from collections import defaultdict

        het_residues = defaultdict(list)

        with open(pdb_file) as f:
            for line in f:
                if line.startswith("HETATM"):
                    resname = line[17:20].strip()
                    chain, resid = parse_chain_resid(line)
                    key = (resname, chain, resid)
                    if resname not in standard_residues:
                        het_residues[key].append(line)

        if not het_residues:
            print("No ligand found in PDB file")
            return None, None

        # Find largest hetero group
        largest = max(het_residues.items(), key=lambda x: len(x[1]))
        (resname, chain, resid), ligand_atoms = largest
        ligand_resname = resname
    else:
        # Extract specified ligand
        ligand_atoms = []
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("HETATM") and line[17:20].strip() == ligand_resname:
                    ligand_atoms.append(line)

        if not ligand_atoms:
            print(f"Ligand {ligand_resname} not found in PDB file")
            return None, None

    # Write ligand PDB
    ligand_pdb = f"{output_prefix}_ligand.pdb"

    # Get atom serial numbers for CONECT records
    ligand_serials = set()
    for line in ligand_atoms:
        try:
            serial = int(line[6:11])
            ligand_serials.add(serial)
        except Exception:
            continue

    # Write ligand PDB with relevant CONECT records
    with open(pdb_file) as inp, open(ligand_pdb, "w") as out:
        # Copy header lines
        for line in inp:
            if line.startswith(("HEADER", "TITLE", "REMARK")):
                out.write(line)
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                break

        # Write ligand atoms
        for line in ligand_atoms:
            out.write(line)

        # Write relevant CONECT records
        inp.seek(0)
        for line in inp:
            if line.startswith("CONECT"):
                fields = line.split()
                if len(fields) < 2:
                    continue
                try:
                    main_serial = int(fields[1])
                except Exception:
                    continue
                if main_serial in ligand_serials:
                    filtered_fields = [fields[0], fields[1]] + [
                        f for f in fields[2:] if int(f) in ligand_serials
                    ]
                    if len(filtered_fields) > 2:
                        out.write("{:<6}".format(filtered_fields[0]))
                        for f in filtered_fields[1:]:
                            out.write("{:>5}".format(f))
                        out.write("\n")

        out.write("END\n")

    print(f"Extracted ligand {ligand_resname} to {ligand_pdb}")
    return ligand_resname, ligand_pdb


def prepare_dna(pdb_file, output_prefix="output"):
    """Extract and prepare DNA structure"""
    dna_residues = ["DA", "DT", "DG", "DC", "A", "T", "G", "C"]
    dna_pdb = f"{output_prefix}_dna_only.pdb"

    # Extract DNA residues
    with open(pdb_file) as inp, open(dna_pdb, "w") as out:
        for line in inp:
            if line.startswith("ATOM"):
                out.write(line)
            elif line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname in dna_residues:
                    out.write(line)
            else:
                out.write(line)

    # Fix missing atoms and add hydrogens
    fixer = PDBFixer(filename=dna_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    fixed_pdb = f"{output_prefix}_dna_fixed.pdb"
    with open(fixed_pdb, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"Prepared DNA structure: {fixed_pdb}")
    return fixed_pdb


def fix_ligand_hydrogens(ligand_pdb, output_pdb):
    """Add missing hydrogens to ligand"""
    fixer = PDBFixer(filename=ligand_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    with open(output_pdb, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"Added hydrogens to ligand: {output_pdb}")
    return output_pdb


def count_electrons(ligand_pdb):
    """Count electrons in ligand to check for odd electron count"""
    electron_count = 0

    with open(ligand_pdb) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            # Try to get element from element column
            elem = line[76:78].strip().upper()
            if not elem:
                # Try to get from atom name
                atom_name = line[12:16].strip().upper()
                if atom_name[:2] in PERIODIC_TABLE:
                    elem = atom_name[:2]
                elif atom_name[:1] in PERIODIC_TABLE:
                    elem = atom_name[:1]

            if elem in PERIODIC_TABLE:
                electron_count += PERIODIC_TABLE[elem]

    return electron_count


def parameterize_ligand(ligand_pdb, ligand_resname, charge=0):
    """Parameterize ligand using AmberTools"""
    # Check electron count
    n_electrons = count_electrons(ligand_pdb)
    n_electrons_total = n_electrons - charge

    if n_electrons_total % 2 != 0:
        raise ValueError(
            f"Ligand {ligand_resname} has an odd number of electrons ({n_electrons_total}). "
            f"Please check the ligand charge (currently {charge})."
        )

    # Run antechamber
    mol2_file = f"{ligand_resname}.mol2"
    frcmod_file = f"{ligand_resname}.frcmod"

    run_command(
        f"antechamber -i {ligand_pdb} -fi pdb -o {mol2_file} -fo mol2 "
        f"-c bcc -nc {charge} -s 2 -at gaff2 -rn {ligand_resname}"
    )

    run_command(f"parmchk2 -i {mol2_file} -f mol2 -o {frcmod_file}")

    print(f"Parameterized ligand: {mol2_file}, {frcmod_file}")
    return mol2_file, frcmod_file


def build_amber_system(
    dna_pdb, ligand_mol2, ligand_frcmod, output_prefix, box_padding=1.0
):
    """Build and solvate system using tleap"""
    buffer_angstrom = float(box_padding) * 10.0  # Convert nm to Å

    leap_script = f"""
source leaprc.gaff2
source leaprc.DNA.OL15
source leaprc.water.tip3p

set default FlexibleWater off

LIG = loadmol2 {ligand_mol2}
loadamberparams {ligand_frcmod}
REC = loadpdb {dna_pdb}
complex = combine {{REC LIG}}
solvateBox complex TIP3PBOX {buffer_angstrom:.2f}
addions complex Na+ 0
saveamberparm complex {output_prefix}.prmtop {output_prefix}.inpcrd
savepdb complex {output_prefix}_leap.pdb
quit
"""

    leap_file = f"{output_prefix}_leap.in"
    with open(leap_file, "w") as f:
        f.write(leap_script)

    run_command(f"tleap -f {leap_file}")

    print(f"Built solvated system: {output_prefix}.prmtop, {output_prefix}.inpcrd")
    return f"{output_prefix}.prmtop", f"{output_prefix}.inpcrd"


def create_openmm_system(
    prmtop_file, inpcrd_file, output_prefix, ml_potential=False, ml_model_path=None
):
    """Create OpenMM system from Amber files"""
    # Load Amber files
    prmtop = app.AmberPrmtopFile(prmtop_file)
    inpcrd = app.AmberInpcrdFile(inpcrd_file)

    # Create system
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )

    # Optionally replace with ML potential
    if ml_potential:
        try:
            from torchmdnet.openmm import TorchMDForce
        except ImportError:
            raise ImportError(
                "TorchMD-NET is not installed. Please install it to use ML potentials."
            )

        if not ml_model_path:
            raise ValueError("ML model path must be provided when using ML potentials.")

        # Remove classical forces
        forces_to_remove = []
        for i, force in enumerate(system.getForces()):
            if force.__class__.__name__ in [
                "HarmonicBondForce",
                "HarmonicAngleForce",
                "PeriodicTorsionForce",
                "NonbondedForce",
            ]:
                forces_to_remove.append(i)

        for idx in sorted(forces_to_remove, reverse=True):
            system.removeForce(idx)

        # Add ML force
        torchmd_force = TorchMDForce(ml_model_path)
        system.addForce(torchmd_force)
        print(f"Added ML potential from {ml_model_path}")

    # Save system files
    with open(f"{output_prefix}_system.xml", "w") as f:
        f.write(mm.XmlSerializer.serialize(system))

    # Create and save integrator
    integrator = mm.LangevinIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    with open(f"{output_prefix}_integrator.xml", "w") as f:
        f.write(mm.XmlSerializer.serialize(integrator))

    # Save structure PDB
    modeller = app.Modeller(prmtop.topology, inpcrd.positions)
    with open(f"{output_prefix}_structure.pdb", "w") as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

    print(f"Created OpenMM system: {output_prefix}_system.xml")
    return f"{output_prefix}_system.xml", f"{output_prefix}_structure.pdb"


def build_system(
    pdb_id,
    output_prefix="output",
    box_padding=1.0,
    ml_potential=False,
    ml_model_path=None,
    ligand_resname=None,
    ligand_charge=0,
    pdb_is_file=False,
):
    """Main function to build a solvated DNA/ligand system"""
    # Get PDB file
    if pdb_is_file:
        pdb_file = pdb_id
    else:
        pdb_file = download_pdb(pdb_id)

    # Extract and prepare ligand
    ligand_resname, ligand_pdb = extract_ligand(pdb_file, ligand_resname, output_prefix)
    if not ligand_pdb:
        raise ValueError("Failed to extract ligand from PDB file")

    # Fix ligand hydrogens
    ligand_pdb_h = f"{output_prefix}_ligand_h.pdb"
    ligand_pdb_h = fix_ligand_hydrogens(ligand_pdb, ligand_pdb_h)

    # Prepare DNA structure
    dna_pdb = prepare_dna(pdb_file, output_prefix)

    # Parameterize ligand
    mol2_file, frcmod_file = parameterize_ligand(
        ligand_pdb_h, ligand_resname, ligand_charge
    )

    # Build Amber system
    prmtop_file, inpcrd_file = build_amber_system(
        dna_pdb, mol2_file, frcmod_file, output_prefix, box_padding
    )

    # Create OpenMM system
    system_xml, structure_pdb = create_openmm_system(
        prmtop_file, inpcrd_file, output_prefix, ml_potential, ml_model_path
    )

    print(f"System preparation complete: {structure_pdb}, {system_xml}")
    return structure_pdb, system_xml
