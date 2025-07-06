# define_binding_site.py

import argparse
import json
from Bio.PDB import PDBParser, NeighborSearch, Selection
import numpy as np

def get_binding_site(pdb_file, ligand_resname, cutoff=5.0, output_json="binding_site.json"):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)

    # Collect all atoms
    atoms = list(structure.get_atoms())

    # Identify ligand atoms
    ligand_atoms = [
        atom for atom in atoms
        if atom.get_parent().get_resname() == ligand_resname
        and atom.get_parent().id[0] == "H_"  # Heteroatom indicator for ligand
    ]
    if not ligand_atoms:
        print(f"No ligand atoms found with residue name {ligand_resname}.")
        return

    # Build NeighborSearch object
    neighbor_search = NeighborSearch(atoms)

    # Find all residues within cutoff of any ligand atom
    nearby_residues = set()
    for ligand_atom in ligand_atoms:
        neighbors = neighbor_search.search(ligand_atom.coord, cutoff, level='R')
        for res in neighbors:
            # Skip the ligand itself
            if res.get_resname() == ligand_resname and res.id[0] == "H_":
                continue
            # (chain, resseq, icode)
            residue_id = (
                res.get_parent().id,  # Chain ID
                res.get_id()[1],      # Residue number
                res.get_resname()
            )
            nearby_residues.add(residue_id)

    # Output as sorted list
    binding_site = [
        {"chain": chain, "resid": resid, "resname": resname}
        for (chain, resid, resname) in sorted(nearby_residues)
    ]

    # Print and save
    print(f"Residues within {cutoff} Å of {ligand_resname}:")
    for res in binding_site:
        print(f"  Chain {res['chain']} Residue {res['resid']} ({res['resname']})")
    with open(output_json, "w") as f:
        json.dump(binding_site, f, indent=2)
    print(f"Binding site saved to {output_json}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Define binding site residues near a ligand.")
    parser.add_argument("--pdb", required=True, help="Input PDB file (from build.py)")
    parser.add_argument("--ligand", required=True, help="Ligand residue name (e.g., DAN)")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Distance cutoff in Å")
    parser.add_argument("--out", default="binding_site.json", help="Output JSON file")
    args = parser.parse_args()
    get_binding_site(args.pdb, args.ligand, args.cutoff, args.out)

