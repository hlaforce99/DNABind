#!/usr/bin/env python3
"""
utils.py

Utility functions for the DNABind pipeline:
- File and directory operations
- Atom selection and filtering
- Unit handling
- Data conversion
- Dependency checking
"""

import os
import json
import glob
import shutil
import subprocess
import re
import numpy as np
from Bio.PDB import PDBParser

# Physical constants
kB = 0.0083144621  # Boltzmann constant in kJ/mol/K


def create_directory(directory_path):
    """Create directory if it doesn't exist"""
    if directory_path and directory_path != "." and not os.path.exists(directory_path):
        os.makedirs(directory_path, exist_ok=True)
        print(f"Created directory: {directory_path}")
    return directory_path


def load_binding_site(binding_site_file):
    """Load binding site definition from JSON file"""
    try:
        with open(binding_site_file, "r") as f:
            binding_site = json.load(f)
        return binding_site
    except Exception as e:
        print(f"Error loading binding site file {binding_site_file}: {e}")
        return []


def get_binding_site_indices(pdb_file, binding_site):
    """
    Get atom indices for binding site residues using Biopython.
    Args:
        pdb_file: Path to PDB file
        binding_site: List of dicts with chain, resid, (optional) resname
    Returns:
        List of atom indices
    """
    if not pdb_file or not os.path.exists(pdb_file):
        print("Error: PDB file required for binding site detection")
        return []

    indices = []
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("complex", pdb_file)
        atom_map = {}
        atom_index = 0
        for chain in structure.get_chains():
            chain_id = chain.id
            for residue in chain:
                if residue.id[0] == " ":
                    resid = residue.id[1]
                    resname = residue.resname
                    key1 = (chain_id, resid, resname)
                    key2 = (chain_id, resid)
                    key3 = (str(chain_id), resid)
                    key4 = (chain_id, str(resid))
                    res_atoms = []
                    for atom in residue:
                        res_atoms.append(atom_index)
                        atom_index += 1
                    atom_map[key1] = res_atoms
                    atom_map[key2] = res_atoms
                    atom_map[key3] = res_atoms
                    atom_map[key4] = res_atoms
                else:
                    for atom in residue:
                        atom_index += 1
        print(f"Looking for binding site residues in PDB...")
        for residue in binding_site:
            chain_id = residue["chain"]
            resid = residue["resid"]
            resname = residue.get("resname", None)
            try:
                resid_int = int(resid)
            except (ValueError, TypeError):
                resid_int = resid
            found = False
            for key in [
                (chain_id, resid_int, resname),
                (chain_id, resid_int),
                (str(chain_id), resid_int),
                (chain_id, str(resid_int)),
            ]:
                if key in atom_map:
                    indices.extend(atom_map[key])
                    found = True
                    print(
                        f"  Found matching residue: Chain {chain_id}, ID {resid}, Name {resname}"
                    )
                    break
            if not found:
                print(
                    f"  Warning: Could not find residue {chain_id}:{resid} ({resname}) in PDB"
                )
    except Exception as e:
        print(f"Error processing PDB file with Biopython: {e}")
        return []
    print(f"Found {len(indices)} atoms in binding site")
    return indices


def get_ligand_indices(pdb_file, ligand_resname):
    """
    Get atom indices for ligand residue using Biopython.
    Args:
        pdb_file: Path to PDB file
        ligand_resname: Residue name of ligand
    Returns:
        List of atom indices
    """
    indices = []
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("complex", pdb_file)
        atom_index = 0
        for chain in structure.get_chains():
            for residue in chain:
                for atom in residue:
                    if residue.resname.strip().upper() == ligand_resname.upper():
                        indices.append(atom_index)
                    atom_index += 1
    except Exception as e:
        print(f"Error processing PDB file with Biopython: {e}")
        return []
    print(f"Found {len(indices)} atoms in ligand {ligand_resname}")
    return indices


def filter_atoms_near_centroid(coords, indices, cutoff):
    """
    Filter atoms to keep only those within cutoff of centroid.
    If cutoff <= 0, return all indices unchanged.

    Args:
        coords: Coordinates array [n_atoms, 3]
        indices: List of atom indices to filter
        cutoff: Distance cutoff in nm

    Returns:
        Filtered list of atom indices
    """
    if cutoff <= 0:
        return indices

    # Calculate centroid
    selected_coords = coords[indices]
    centroid = np.mean(selected_coords, axis=0)

    # Calculate distances to centroid
    dists = np.sqrt(np.sum((selected_coords - centroid) ** 2, axis=1))

    # Filter atoms
    filtered_indices = [indices[i] for i in range(len(indices)) if dists[i] <= cutoff]

    return filtered_indices


def check_dependencies():
    """
    Check if required dependencies are installed.

    Returns:
        bool: True if all dependencies are available
    """
    missing = []

    # Try importing required packages
    try:
        import openmm
    except ImportError:
        missing.append("OpenMM")

    try:
        import mdtraj
    except ImportError:
        missing.append("MDTraj")

    try:
        import pdbfixer
    except ImportError:
        missing.append("PDBFixer")

    try:
        import matplotlib
    except ImportError:
        missing.append("Matplotlib")

    try:
        import scipy
    except ImportError:
        missing.append("SciPy")

    try:
        import Bio
    except ImportError:
        missing.append("BioPython")

    # Check for command-line tools
    def check_command(cmd):
        try:
            result = subprocess.run(
                f"which {cmd}",
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            return result.returncode == 0
        except Exception:
            return False

    for cmd in ["antechamber", "parmchk2", "tleap"]:
        if not check_command(cmd):
            missing.append(cmd)

    if missing:
        print("Missing dependencies:")
        for dep in missing:
            print(f"  - {dep}")
        print("\nPlease install missing dependencies before running DNABind.")
        return False

    return True


def to_serializable(obj):
    """Recursively convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [to_serializable(v) for v in obj]
    else:
        return obj


def parse_chain_resid(line):
    """
    Robustly parse chain ID and residue number from a PDB ATOM/HETATM line.
    Handles both properly and improperly spaced formats.
    """
    # PDB format: columns 22 (chain), 23-26 (resid)
    chain = line[21]
    resid_field = line[22:26]
    # If resid_field contains both chain and resid (e.g. 'A1000'), split it
    match = re.match(r"([A-Za-z])?\s*(\d+)", resid_field)
    if match:
        if match.group(1):
            chain = match.group(1)
        resid = int(match.group(2))
    else:
        # Fallback: try to parse as int
        try:
            resid = int(resid_field.strip())
        except Exception:
            resid = None
    return chain, resid


def cleanup_files(output_prefix, output_dir="output", root_dir=".", pdb_id=None):
    """
    Move all important files (keep_files + pdb_id) to the output directory,
    and remove temporary files matching temp_patterns from both the output and root directories.
    """
    import os, shutil, glob

    # --- 1. Define keep_files ---
    keep_files = {
        f"{output_prefix}_equilibrated.pdb",
        f"{output_prefix}_equilibrated_system.xml",
        f"{output_prefix}_equilibrated_integrator.xml",
        f"{output_prefix}_structure.pdb",
        f"{output_prefix}_system.xml",
        f"{output_prefix}_integrator.xml",
        f"{output_prefix}_pmf.pdf",
        f"{output_prefix}_prmtop",
        f"{output_prefix}.inpcrd",
        "binding_site.json",
        "traj.dcd",
        f"{output_prefix}_bias.npz",
        f"{output_prefix}_binding_energy.json",
        f"{output_prefix}_bound_unbound.npy",
    }
    if pdb_id:
        keep_files.add(f"{pdb_id}.pdb")

    # --- 2. Move keep_files to output_dir if needed ---
    for fname in keep_files:
        if not fname:
            continue
        src_path = os.path.abspath(os.path.join(root_dir, fname))
        dst_path = os.path.abspath(os.path.join(output_dir, fname))
        # Only move if src_path exists and is not the same as dst_path
        if src_path != dst_path and os.path.isfile(src_path):
            try:
                shutil.move(src_path, dst_path)
                print(f"Moved {fname} to {output_dir}")
            except Exception as e:
                print(f"Could not move {fname}: {e}")

    # --- 3. Compute absolute paths of keep_files (for later) ---
    keep_files_abs = set()
    for f in keep_files:
        if not f:
            continue
        keep_files_abs.add(os.path.abspath(os.path.join(output_dir, f)))
        keep_files_abs.add(os.path.abspath(os.path.join(root_dir, f)))

    # --- 4. Remove temp files not in keep_files ---
    temp_patterns = [
        "*.mol2",
        "*.frcmod",
        "*.log",
        "ANTECHAMBER*",
        "sqm*",
        f"{output_prefix}_ligand*.pdb",
        f"{output_prefix}_dna*.pdb",
        f"{output_prefix}_leap.*",
        f"{output_prefix}_fixed.pdb",
        "ATOMTYPE.INF",
    ]

    for directory in [output_dir, root_dir]:
        for pattern in temp_patterns:
            for file_path in glob.glob(os.path.join(directory, pattern)):
                abs_path = os.path.abspath(file_path)
                if os.path.isfile(abs_path) and abs_path not in keep_files_abs:
                    try:
                        os.remove(abs_path)
                        print(f"Removed temporary file: {file_path}")
                    except Exception as e:
                        print(f"Could not remove {file_path}: {e}")
