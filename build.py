# build.py

import argparse
import os
import sys
import subprocess
import warnings
from pdbfixer import PDBFixer
import openmm.app as app
import openmm as mm
import openmm.unit as unit

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

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


def atom_element_from_pdb_line(line):
    elem = line[76:78].strip().upper()
    if elem:
        return elem
    name = line[12:16].strip().upper()
    if name[:2] in PERIODIC_TABLE:
        return name[:2]
    elif name[:1] in PERIODIC_TABLE:
        return name[:1]
    return None


def count_ligand_electrons(ligand_atoms):
    n_e = 0
    for line in ligand_atoms:
        elem = atom_element_from_pdb_line(line)
        if elem in PERIODIC_TABLE:
            n_e += PERIODIC_TABLE[elem]
    return n_e


def extract_ligand_atoms(pdb_filename, ligand_resname=None):
    standard = set(["DA", "DT", "DG", "DC", "A", "T", "G", "C", "HOH", "WAT"])
    ligand_lines = []
    if ligand_resname is None:
        from collections import defaultdict

        het_residues = defaultdict(list)
        with open(pdb_filename) as f:
            for line in f:
                if line.startswith("HETATM"):
                    resname = line[17:20].strip()
                    chain = line[21]
                    resid = line[22:26].strip()
                    key = (resname, chain, resid)
                    if resname not in standard:
                        het_residues[key].append(line)
        if not het_residues:
            return None, []
        largest = max(het_residues.items(), key=lambda x: len(x[1]))
        (resname, chain, resid), atoms = largest
        return resname, atoms
    else:
        with open(pdb_filename, "r") as f:
            for line in f:
                if line.startswith("HETATM") and line[17:20].strip() == ligand_resname:
                    ligand_lines.append(line)
        return ligand_resname, ligand_lines


def write_ligand_pdb(ligand_atoms, output_pdb, original_pdb):
    ligand_serials = set()
    for line in ligand_atoms:
        if line.startswith(("HETATM", "ATOM")):
            try:
                serial = int(line[6:11])
                ligand_serials.add(serial)
            except Exception:
                continue
    with open(original_pdb, "r") as inp, open(output_pdb, "w") as out:
        for line in inp:
            if line.startswith(("HEADER", "TITLE", "REMARK")):
                out.write(line)
        for line in ligand_atoms:
            out.write(line)
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
                if main_serial not in ligand_serials:
                    continue
                filtered_fields = [fields[0], fields[1]] + [
                    f for f in fields[2:] if int(f) in ligand_serials
                ]
                if len(filtered_fields) > 2:
                    out.write("{:<6}".format(filtered_fields[0]))
                    for f in filtered_fields[1:]:
                        out.write("{:>5}".format(f))
                    out.write("\n")
        out.write("END\n")


def remove_heterogens(input_pdb, output_pdb, keep_resnames=None):
    if keep_resnames is None:
        keep_resnames = ["DA", "DT", "DG", "DC", "A", "T", "G", "C"]
    with open(input_pdb) as inp, open(output_pdb, "w") as out:
        for line in inp:
            if line.startswith("ATOM"):
                out.write(line)
            elif line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname in keep_resnames:
                    out.write(line)
            else:
                out.write(line)


def run(cmd, **kwargs):
    result = subprocess.run(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kwargs
    )
    if result.returncode != 0:
        sys.stderr.write(result.stdout.decode())
        raise RuntimeError(f"Command failed: {cmd}")
    return result


def fix_ligand_hydrogens(ligand_pdb_in, ligand_pdb_out):
    fixer = PDBFixer(filename=ligand_pdb_in)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(ligand_pdb_out, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)


def parameterize_ligand_with_ambertools(
    ligand_pdb, ligand_resname, workdir=".", charge=0, gaff_version="gaff2"
):
    base = os.path.join(workdir, ligand_resname)
    mol2 = f"{base}.mol2"
    frcmod = f"{base}.frcmod"
    run(
        f"antechamber -i {ligand_pdb} -fi pdb -o {mol2} -fo mol2 -c bcc -nc {charge} -s 2 -at {gaff_version} -rn {ligand_resname}"
    )
    run(f"parmchk2 -i {mol2} -f mol2 -o {frcmod}")
    return mol2, frcmod


def build_leap_input(
    receptor_pdb,
    ligand_mol2,
    ligand_frcmod,
    output_prefix,
    box_padding=1.0,
):
    buffer_angstrom = float(box_padding) * 10.0
    leap = f"""
source leaprc.gaff2
source leaprc.DNA.OL15
source leaprc.water.tip3p

set default FlexibleWater off

LIG = loadmol2 {ligand_mol2}
loadamberparams {ligand_frcmod}
REC = loadpdb {receptor_pdb}
complex = combine {{REC LIG}}
solvateBox complex TIP3PBOX {buffer_angstrom:.2f}
addions complex Na+ 0
saveamberparm complex {output_prefix}.prmtop {output_prefix}.inpcrd
savepdb complex {output_prefix}_leap.pdb
quit
"""
    return leap


def cleanup_temp_files(prefix, ligand_resname):
    exts = [
        ".frcmod",
        ".mol2",
        ".log",
        ".out",
        ".in",
        ".pdb",
        ".leap.in",
        "_leap.pdb",
        "_ligand.pdb",
        "_ligand_h.pdb",
        "_dna_only.pdb",
        "_fixed.pdb",
        ".INF",
    ]
    keep = {
        f"{prefix}_structure.pdb",
        f"{prefix}_system.xml",
        f"{prefix}_integrator.xml",
        f"{prefix}.prmtop",
        f"{prefix}.inpcrd",
    }
    for f in os.listdir("."):
        if (
            any(f.endswith(ext) for ext in exts)
            or f.startswith("ANTECHAMBER")
            or f.startswith("sqm")
        ):
            if f not in keep:
                try:
                    os.remove(f)
                except Exception:
                    pass


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
    if pdb_is_file:
        pdb_filename = pdb_id
    else:
        pdb_filename = f"{pdb_id}.pdb"
        if not os.path.exists(pdb_filename):
            from urllib.request import urlopen

            with open(pdb_filename, "w") as f:
                f.write(
                    urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb")
                    .read()
                    .decode()
                )

    ligand_atoms = []
    ligand_pdb = None
    if ligand_resname:
        resname, ligand_atoms = extract_ligand_atoms(pdb_filename, ligand_resname)
        if not ligand_atoms:
            resname, ligand_atoms = extract_ligand_atoms(pdb_filename, None)
            ligand_resname = resname
    else:
        resname, ligand_atoms = extract_ligand_atoms(pdb_filename, None)
        ligand_resname = resname

    ligand_pdb = f"{output_prefix}_ligand.pdb"
    if ligand_atoms:
        write_ligand_pdb(ligand_atoms, ligand_pdb, pdb_filename)
        ligand_pdb_h = f"{output_prefix}_ligand_h.pdb"
        fix_ligand_hydrogens(ligand_pdb, ligand_pdb_h)
        ligand_atoms_h = []
        with open(ligand_pdb_h) as f:
            for line in f:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    ligand_atoms_h.append(line)
        n_e = count_ligand_electrons(ligand_atoms_h)
        n_e_total = n_e - ligand_charge
        if n_e_total % 2 != 0:
            sys.stderr.write(
                f"ERROR: Ligand {ligand_resname} has an odd number of electrons ({n_e_total}).\n"
            )
            sys.exit(1)
        ligand_pdb = ligand_pdb_h
    else:
        ligand_pdb = None

    keep_resnames = ["DA", "DT", "DG", "DC", "A", "T", "G", "C"]
    dna_only_pdb = f"{output_prefix}_dna_only.pdb"
    remove_heterogens(pdb_filename, dna_only_pdb, keep_resnames=keep_resnames)

    fixer = PDBFixer(filename=dna_only_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    fixed_pdb_filename = f"{output_prefix}_fixed.pdb"
    with open(fixed_pdb_filename, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    if ligand_pdb and os.path.getsize(ligand_pdb) > 0:
        mol2, frcmod = parameterize_ligand_with_ambertools(
            ligand_pdb, ligand_resname, workdir=".", charge=ligand_charge
        )
    else:
        mol2 = None
        frcmod = None

    if mol2 and frcmod:
        leap_script = build_leap_input(
            fixed_pdb_filename,
            mol2,
            frcmod,
            output_prefix,
            box_padding=box_padding,
        )
        leapin = f"{output_prefix}_leap.in"
        with open(leapin, "w") as f:
            f.write(leap_script)
        run(f"tleap -f {leapin}")

        prmtop = app.AmberPrmtopFile(f"{output_prefix}.prmtop")
        inpcrd = app.AmberInpcrdFile(f"{output_prefix}.inpcrd")
        modeller = app.Modeller(prmtop.topology, inpcrd.positions)
    else:
        sys.stderr.write("ERROR: Ligand parameterization failed or ligand not found.\n")
        sys.exit(1)

    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )
    if ml_potential:
        try:
            from torchmdnet.openmm import TorchMDForce
        except ImportError:
            sys.stderr.write("ERROR: TorchMD-NET is not installed.\n")
            sys.exit(1)
        if not ml_model_path:
            sys.stderr.write(
                "ERROR: Must provide --ml_model_path when --ml_potential is enabled.\n"
            )
            sys.exit(1)
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
        torchmd_force = TorchMDForce(ml_model_path)
        system.addForce(torchmd_force)

    structure_pdb = f"{output_prefix}_structure.pdb"
    with open(structure_pdb, "w") as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
    with open(f"{output_prefix}_system.xml", "w") as f:
        f.write(mm.XmlSerializer.serialize(system))
    with open(f"{output_prefix}_integrator.xml", "w") as f:
        integrator = mm.LangevinIntegrator(300, 1.0, 0.002)
        f.write(mm.XmlSerializer.serialize(integrator))

    cleanup_temp_files(output_prefix, ligand_resname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build solvated DNA/ligand system for OpenMM using AmberTools for ligand."
    )
    parser.add_argument("--pdb", required=True, help="PDB ID or file")
    parser.add_argument("--prefix", default="output", help="Output prefix")
    parser.add_argument("--box", type=float, default=1.0, help="Box padding (nm)")
    parser.add_argument(
        "--ml_potential", action="store_true", help="Use ML potential (TorchMD-NET)"
    )
    parser.add_argument("--ml_model_path", default=None, help="TorchMD-NET model path")
    parser.add_argument("--ligand", default=None, help="Ligand residue name to extract")
    parser.add_argument(
        "--ligand_charge",
        type=int,
        default=0,
        help="Net charge of ligand (required for correct AM1-BCC charge assignment)",
    )
    parser.add_argument(
        "--pdb_is_file", action="store_true", help="Treat --pdb as a local file"
    )
    args = parser.parse_args()
    build_system(
        pdb_id=args.pdb,
        output_prefix=args.prefix,
        box_padding=args.box,
        ml_potential=args.ml_potential,
        ml_model_path=args.ml_model_path,
        ligand_resname=args.ligand,
        ligand_charge=args.ligand_charge,
        pdb_is_file=args.pdb_is_file,
    )
