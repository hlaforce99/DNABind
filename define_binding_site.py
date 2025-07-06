# define_binding_site.py

import argparse
import json
from Bio.PDB import PDBParser, NeighborSearch


def get_binding_site_auto(pdb_file, ligand_resname, cutoff=5.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    atoms = list(structure.get_atoms())

    ligand_atoms = [
        atom
        for atom in atoms
        if atom.get_parent().get_resname() == ligand_resname
        and atom.get_parent().id[0] == "H_"
    ]
    if not ligand_atoms:
        print(f"No ligand atoms found with residue name {ligand_resname}.")
        return []

    neighbor_search = NeighborSearch(atoms)
    nearby_residues = set()
    for ligand_atom in ligand_atoms:
        neighbors = neighbor_search.search(ligand_atom.coord, cutoff, level="R")
        for res in neighbors:
            if res.get_resname() == ligand_resname and res.id[0] == "H_":
                continue
            residue_id = (
                res.get_parent().id,  # Chain
                res.get_id()[1],  # Residue number
                res.get_resname(),
            )
            nearby_residues.add(residue_id)

    binding_site = [
        {"chain": chain, "resid": resid, "resname": resname}
        for (chain, resid, resname) in sorted(nearby_residues)
    ]
    return binding_site


def get_binding_site_manual(residues):
    # residues: list of "chain:resnum" (e.g., ["A:10", "B:15"])
    binding_site = []
    for r in residues:
        try:
            chain, resid = r.split(":")
            resid = int(resid)
            binding_site.append({"chain": chain, "resid": resid, "resname": None})
        except Exception as e:
            print(f"Could not parse residue {r}: {e}")
    return binding_site


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Define binding site residues.")
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--ligand", help="Ligand residue name (for auto mode)")
    parser.add_argument(
        "--cutoff", type=float, default=5.0, help="Distance cutoff in Å (auto mode)"
    )
    parser.add_argument(
        "--residues",
        nargs="+",
        help="Manual list of residues as chain:resnum (e.g. A:10 B:15)",
    )
    parser.add_argument("--out", default="binding_site.json", help="Output JSON file")
    args = parser.parse_args()

    if args.residues:
        binding_site = get_binding_site_manual(args.residues)
        print("Manual binding site selection:")
    elif args.ligand:
        binding_site = get_binding_site_auto(args.pdb, args.ligand, args.cutoff)
        print(
            f"Automatic binding site selection: residues within {args.cutoff} Å of {args.ligand}"
        )
    else:
        print("Error: must specify either --ligand or --residues")
        exit(1)

    for res in binding_site:
        print(res)
    with open(args.out, "w") as f:
        json.dump(binding_site, f, indent=2)
    print(f"Binding site saved to {args.out}")
