#!/usr/bin/env python3
"""
binding_site.py

Module for defining binding site residues for molecular simulations:
- Can automatically detect residues near a ligand (auto mode)
- Or accept manual residue specification (manual mode)
- Excludes solvent and common ions from the binding site
"""

import json
from Bio.PDB import PDBParser, NeighborSearch

# Common solvent and ion residue names to exclude from binding site
SOLVENT_IONS = {"HOH", "WAT", "NA+", "CL-", "K+", "MG2+", "CA2+", "SO42-", "PO43-"}


def get_binding_site_auto(pdb_file, ligand_resname, cutoff=4.0):
    """
    Automatically detect binding site residues within `cutoff` Å of any ligand atom.
    Excludes the ligand itself and solvent/ions.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    atoms = list(structure.get_atoms())

    # Find ligand atoms
    ligand_atoms = [
        atom for atom in atoms if atom.get_parent().get_resname() == ligand_resname
    ]

    if not ligand_atoms:
        print(f"Warning: Ligand {ligand_resname} not found in {pdb_file}")
        return []

    # Set up neighbor search
    neighbor_search = NeighborSearch(atoms)
    nearby_residues = set()

    # Find residues near ligand atoms
    for ligand_atom in ligand_atoms:
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
    print(
        f"Found {len(binding_site)} residues within {cutoff}Å of ligand {ligand_resname}"
    )
    return binding_site


def get_binding_site_manual(residues, pdb_file=None):
    """
    Accept a manual list of residues in the form chain:resid (e.g. A:10).
    Optionally looks up residue names from PDB file.
    """
    binding_site = []

    # Get residue names from PDB if available
    resname_lookup = {}
    if pdb_file:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", pdb_file)
        for chain in structure.get_chains():
            chain_id = chain.id
            for residue in chain:
                if residue.id[0] == " ":  # Standard residue
                    resid = residue.id[1]
                    resname = residue.resname
                    resname_lookup[(chain_id, resid)] = resname

    # Process manual residue list
    for r in residues:
        try:
            chain, resid = r.split(":")
            resid = int(resid)

            # Look up residue name if available
            resname = resname_lookup.get((chain, resid), None)
            binding_site.append({"chain": chain, "resid": resid, "resname": resname})
        except Exception as e:
            print(f"Warning: Could not parse residue specification '{r}': {e}")
            continue

    return binding_site


def define_binding_site(
    pdb_file,
    ligand_resname=None,
    residues=None,
    cutoff=3.5,
    output_file="binding_site.json",
):
    """
    Define binding site residues, either automatically or manually.

    Args:
        pdb_file: Path to PDB file
        ligand_resname: Residue name of ligand (for auto mode)
        residues: List of residues as "chain:resid" (for manual mode)
        cutoff: Distance cutoff in Å (for auto mode)
        output_file: Output JSON file path

    Returns:
        List of binding site residues as dicts
    """
    # Choose mode: manual or auto
    if residues:
        binding_site = get_binding_site_manual(residues, pdb_file)
        print(f"Using manually specified binding site: {len(binding_site)} residues")
    elif ligand_resname:
        binding_site = get_binding_site_auto(pdb_file, ligand_resname, cutoff)
        print(
            f"Using automatically detected binding site: {len(binding_site)} residues"
        )
    else:
        raise ValueError("Must specify either ligand_resname or residues")

    # Save binding site definition to JSON
    if output_file:
        with open(output_file, "w") as f:
            json.dump(binding_site, f, indent=2)
        print(f"Binding site definition saved to {output_file}")

    return binding_site
