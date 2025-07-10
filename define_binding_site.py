#!/usr/bin/env python3
"""
define_binding_site.py

Utility to define binding site residues for molecular simulations.
- Can automatically detect residues near a ligand (auto mode)
- Or accept manual residue specification (manual mode)
- Excludes solvent and common ions from the binding site

Outputs a JSON file describing the binding site for downstream analysis.
"""

import argparse
import json
from Bio.PDB import PDBParser, NeighborSearch
import sys

# Common solvent and ion residue names to exclude from binding site
SOLVENT_IONS = {"HOH", "WAT", "NA", "CL", "K", "MG", "CA", "SO4", "PO4"}


def get_binding_site_auto(pdb_file, ligand_resname, cutoff=5.0):
    """
    Automatically detect binding site residues within `cutoff` Å of any ligand atom.
    Excludes the ligand itself and solvent/ions.

    Args:
        pdb_file (str): Path to PDB file.
        ligand_resname (str): Residue name of the ligand.
        cutoff (float): Distance cutoff in Ångstroms.

    Returns:
        List[dict]: List of binding site residues as dicts with chain, resid, resname.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    atoms = list(structure.get_atoms())
    ligand_atoms = [
        atom for atom in atoms if atom.get_parent().get_resname() == ligand_resname
    ]
    if not ligand_atoms:
        return []
    neighbor_search = NeighborSearch(atoms)
    nearby_residues = set()
    for ligand_atom in ligand_atoms:
        # Find all residues within cutoff Å of this ligand atom
        neighbors = neighbor_search.search(ligand_atom.coord, cutoff, level="R")
        for res in neighbors:
            # Skip the ligand residue itself
            if res.get_resname() == ligand_resname:
                continue
            residue_id = (
                res.get_parent().id,  # Chain ID
                res.get_id()[1],  # Residue number
                res.get_resname(),  # Residue name
            )
            nearby_residues.add(residue_id)
    # Exclude solvent and ions
    binding_site = [
        {"chain": chain, "resid": resid, "resname": resname}
        for (chain, resid, resname) in sorted(nearby_residues)
        if resname.strip().upper() not in SOLVENT_IONS
    ]
    return binding_site


def get_binding_site_manual(residues):
    """
    Accept a manual list of residues in the form chain:resid (e.g. A:10).
    Residue name is set to None.

    Args:
        residues (List[str]): List like ["A:10", "B:15"]

    Returns:
        List[dict]: List of residue dicts with chain, resid, resname=None.
    """
    binding_site = []
    for r in residues:
        try:
            chain, resid = r.split(":")
            resid = int(resid)
            binding_site.append({"chain": chain, "resid": resid, "resname": None})
        except Exception:
            continue  # Skip malformed residue specs
    return binding_site


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Define binding site residues for molecular simulations."
    )
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--ligand", help="Ligand residue name (for auto mode)")
    parser.add_argument(
        "--cutoff", type=float, default=3.5, help="Distance cutoff in Å (auto mode)"
    )
    parser.add_argument(
        "--residues",
        nargs="+",
        help="Manual list of residues as chain:resnum (e.g. A:10 B:15)",
    )
    parser.add_argument("--out", default="binding_site.json", help="Output JSON file")
    args = parser.parse_args()

    # Choose mode: manual or auto
    if args.residues:
        binding_site = get_binding_site_manual(args.residues)
    elif args.ligand:
        binding_site = get_binding_site_auto(args.pdb, args.ligand, args.cutoff)
    else:
        sys.stderr.write("Error: must specify either --ligand or --residues\n")
        exit(1)

    # Save binding site definition to JSON
    with open(args.out, "w") as f:
        json.dump(binding_site, f, indent=2)
