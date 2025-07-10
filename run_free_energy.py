#!/usr/bin/env python3
"""
run_free_energy.py

A pipeline to perform ligand-DNA binding/unbinding free energy calculations using OpenMM and metadynamics.
- Sets up the system, runs energy minimization and equilibration,
- Applies wall and COM restraints,
- Runs production metadynamics MD,
- Analyzes the resulting trajectory to estimate binding free energy (ΔG).

"""

import argparse
import json
import numpy as np
import time
import sys

from openmm.app import PDBFile, Simulation, StateDataReporter, DCDReporter
from openmm import (
    XmlSerializer,
    unit,
    CustomCentroidBondForce,
    LangevinIntegrator,
    MonteCarloBarostat,
    CustomExternalForce,
)
from openmm.unit import kelvin, bar
import mdtraj as md
from openmm.app.metadynamics import Metadynamics, BiasVariable


def strip_unit(x):
    """
    Convert OpenMM Quantity objects to numpy arrays in nanometers.
    Handles both Quantity and plain arrays.
    """
    try:
        return np.array(x.value_in_unit(unit.nanometer))
    except AttributeError:
        return np.array(x)


def load_binding_site(json_file):
    """
    Load binding site residue definitions from a JSON file.
    """
    with open(json_file, "r") as f:
        return json.load(f)


def is_heavy_atom(atom):
    """
    Returns True if the atom is not a hydrogen (heavy atom).
    """
    return atom.element.symbol != "H"


def get_binding_site_indices(topology, binding_site):
    """
    Given an OpenMM topology and a list of binding site residues (dicts),
    return the indices of all heavy atoms in the site.
    """
    md_top = md.Topology.from_openmm(topology)
    table, _ = md_top.to_dataframe()
    chain_map = {c.chain_id: i for i, c in enumerate(md_top.chains)}
    indices = []
    for res in binding_site:
        chain_index = chain_map.get(res["chain"])
        if chain_index is None:
            continue
        matches = table[
            (table["chainID"] == chain_index)
            & (table["resSeq"] == int(res["resid"]))
            & (table["resName"] == res["resname"])
        ]
        for idx in matches.index.tolist():
            atom = list(md_top.atoms)[idx]
            if is_heavy_atom(atom):
                indices.append(idx)
    return indices


def get_ligand_indices(topology, ligand_resname):
    """
    Return all heavy atom indices for the ligand, given residue name.
    """
    atoms = list(topology.atoms()) if callable(topology.atoms) else topology.atoms
    return [
        atom.index
        for atom in atoms
        if atom.residue.name == ligand_resname and is_heavy_atom(atom)
    ]


def get_centroid(positions, indices):
    """
    Compute the centroid of a group of atoms.
    """
    coords = np.array([strip_unit(positions[i]) for i in indices])
    return coords.mean(axis=0)


def filter_atoms_near_centroid(positions, indices, cutoff_nm):
    """
    Select only those atoms within cutoff_nm of the centroid of the group.
    Used to focus the CV on the 'core' of each group.
    """
    if not indices:
        return []
    centroid = get_centroid(positions, indices)
    coords = np.array([strip_unit(positions[i]) for i in indices])
    dists = np.linalg.norm(coords - centroid, axis=1)
    return [idx for idx, d in zip(indices, dists) if d <= cutoff_nm]


def create_centroid_distance_cv(ligand_indices, site_indices):
    """
    Create a CustomCentroidBondForce to measure the centroid-centroid distance
    between ligand and binding site.
    """
    force = CustomCentroidBondForce(2, "distance(g1,g2)")
    force.addGroup(ligand_indices)
    force.addGroup(site_indices)
    force.addBond([0, 1], [])
    return force


def create_ligand_wall_restraint(
    ligand_indices, box_vectors, wall_buffer=2.0, k_wall=100.0
):
    """
    Returns a CustomCentroidBondForce that restrains the centroid of a group
    (e.g., ligand) to stay inside the solvent box, with a soft wall.
    """
    box_lengths = [box_vectors[i][i]._value for i in range(3)]
    lower = [wall_buffer] * 3
    upper = [box_lengths[i] - wall_buffer for i in range(3)]
    expr = (
        "step(xmin-x1)*k_wall*(xmin-x1)^2 + step(x1-xmax)*k_wall*(x1-xmax)^2 + "
        "step(ymin-y1)*k_wall*(ymin-y1)^2 + step(y1-ymax)*k_wall*(y1-ymax)^2 + "
        "step(zmin-z1)*k_wall*(zmin-z1)^2 + step(z1-zmax)*k_wall*(z1-zmax)^2"
    )
    force = CustomCentroidBondForce(1, expr)
    force.addGroup(ligand_indices)
    force.addPerBondParameter("xmin")
    force.addPerBondParameter("xmax")
    force.addPerBondParameter("ymin")
    force.addPerBondParameter("ymax")
    force.addPerBondParameter("zmin")
    force.addPerBondParameter("zmax")
    force.addPerBondParameter("k_wall")
    force.addBond([0], [*lower, *upper, k_wall])
    return force


def add_com_restraint(system, atom_indices, k_rest=0.01):
    """
    Add a weak harmonic restraint to keep the center of mass of a group
    (e.g., DNA) near the box center.
    """

    expr = "k_rest*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
    force = CustomExternalForce(expr)
    box_center = [
        0.5 * system.getDefaultPeriodicBoxVectors()[i][i]._value for i in range(3)
    ]
    for idx in atom_indices:
        force.addParticle(idx, [k_rest, *box_center])
    system.addForce(force)


def run_md_with_metadynamics(
    system_xml,
    pdb_file,
    ligand_indices,
    site_indices,
    nsteps=500000,
    report_interval=1000,
    output_traj="traj.dcd",
    metad_height=10.0,
    metad_sigma_dist=0.2,
    metad_biasfactor=10,
    metad_stride=1000,
    temperature=300,
    pressure=1.0,
    wall_buffer=2.0,
    k_wall=100.0,
    nvt_steps=10000,
    npt_steps=20000,
    eq_steps=10000,
    output_prefix="output",
):
    """
    Main function for running the MD protocol:
    - Minimization
    - NVT and NPT equilibration
    - Metadynamics production
    - Applies wall and COM restraints
    Returns: (trajectory_file, topology, equilibrated_pdb_filename)
    """
    # Load system and structure
    system = XmlSerializer.deserialize(open(system_xml).read())
    pdb = PDBFile(pdb_file)
    box_vectors = pdb.topology.getPeriodicBoxVectors()

    # Wall restraints for ligand and DNA
    wall_force_ligand = create_ligand_wall_restraint(
        ligand_indices, box_vectors, wall_buffer, k_wall
    )
    dna_resnames = {"DA", "DT", "DG", "DC"}
    dna_indices = [
        atom.index for atom in pdb.topology.atoms() if atom.residue.name in dna_resnames
    ]
    wall_force_dna = create_ligand_wall_restraint(
        dna_indices, box_vectors, wall_buffer, k_wall / 5.0
    )
    add_com_restraint(system, dna_indices, k_rest=0.01)
    system.addForce(wall_force_ligand)
    system.addForce(wall_force_dna)

    # --- Minimization ---
    integrator = LangevinIntegrator(
        temperature * kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    print("Energy minimization...")
    simulation.minimizeEnergy(maxIterations=5000)

    # --- NVT Equilibration ---
    simulation.context.setVelocitiesToTemperature(temperature * kelvin)
    print(f"NVT equilibration ({nvt_steps} steps)...")
    simulation.step(nvt_steps)

    # --- NPT Equilibration ---
    barostat = MonteCarloBarostat(pressure * bar, temperature * kelvin)
    system.addForce(barostat)
    integrator_npt = LangevinIntegrator(
        temperature * kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation_npt = Simulation(pdb.topology, system, integrator_npt)
    # Carry over positions/velocities with enforcePeriodicBox=True!
    simulation_npt.context.setPositions(
        simulation.context.getState(
            getPositions=True, enforcePeriodicBox=True
        ).getPositions()
    )
    simulation_npt.context.setVelocities(
        simulation.context.getState(getVelocities=True).getVelocities()
    )
    print(f"NPT equilibration ({npt_steps} steps)...")
    simulation_npt.reporters.append(
        StateDataReporter(
            sys.stdout,
            1000,
            step=True,
            potentialEnergy=True,
            temperature=True,
            volume=True,
            density=True,
        )
    )
    simulation_npt.step(npt_steps)
    # Remove barostat for production
    for i in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(i), MonteCarloBarostat):
            system.removeForce(i)
            break

    # --- Final NVT Equilibration ---
    state = simulation_npt.context.getState(
        getPositions=True, getVelocities=True, enforcePeriodicBox=True
    )
    positions = state.getPositions()
    velocities = state.getVelocities()
    box_vectors = state.getPeriodicBoxVectors()
    integrator_eq = LangevinIntegrator(
        temperature * kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation_eq = Simulation(pdb.topology, system, integrator_eq)
    simulation_eq.context.setPositions(positions)
    simulation_eq.context.setVelocities(velocities)
    simulation_eq.context.setPeriodicBoxVectors(*box_vectors)
    print(f"Final NVT equilibration ({eq_steps} steps)...")
    simulation_eq.step(eq_steps)

    # Save equilibrated structure to PDB file
    eq_state = simulation_eq.context.getState(
        getPositions=True, enforcePeriodicBox=True
    )
    eq_positions = eq_state.getPositions()
    eq_box = eq_state.getPeriodicBoxVectors()
    equilibrated_pdb = f"{output_prefix}_equilibrated.pdb"
    pdb.topology.setPeriodicBoxVectors(eq_box)
    with open(equilibrated_pdb, "w") as f:
        PDBFile.writeFile(pdb.topology, eq_positions, f, keepIds=True)
    print(f"Saved equilibrated structure to {equilibrated_pdb}")

    # --- Metadynamics Setup ---
    # The CV is the centroid distance between ligand and binding site
    centroid_dist_force = create_centroid_distance_cv(ligand_indices, site_indices)
    centroid_dist_var = BiasVariable(
        centroid_dist_force,
        minValue=0.0,
        maxValue=5.0,
        biasWidth=metad_sigma_dist,
        periodic=False,
        gridWidth=200,
    )
    meta = Metadynamics(
        system,
        [centroid_dist_var],
        temperature * unit.kelvin,
        biasFactor=metad_biasfactor,
        height=metad_height * unit.kilojoule_per_mole,
        frequency=metad_stride,
        saveFrequency=report_interval,
        biasDir=".",
    )

    # --- Production Simulation with Metadynamics ---
    integrator_prod = LangevinIntegrator(
        temperature * kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation_prod = Simulation(pdb.topology, system, integrator_prod)
    simulation_prod.context.setPositions(
        simulation_eq.context.getState(
            getPositions=True, enforcePeriodicBox=True
        ).getPositions()
    )
    simulation_prod.context.setVelocities(
        simulation_eq.context.getState(getVelocities=True).getVelocities()
    )
    simulation_prod.context.setPeriodicBoxVectors(*box_vectors)
    dcd_stride = max(1000, report_interval)
    simulation_prod.reporters.append(DCDReporter(output_traj, dcd_stride))
    simulation_prod.reporters.append(
        StateDataReporter(
            sys.stdout,
            report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=nsteps,
            separator="\t",
        )
    )
    print(
        f"Production (NVT+MetaD): {nsteps} steps, hill height {metad_height}, sigma {metad_sigma_dist}, stride {metad_stride}, bias factor {metad_biasfactor}"
    )
    meta.step(simulation_prod, nsteps)
    return output_traj, pdb.topology, equilibrated_pdb


def block_average_deltaG(bound, kT=2.479, n_blocks=5):
    """
    Estimate ΔG and error bars via block averaging of the bound/unbound state array.
    Returns means and standard deviations for ΔG, P_bound, P_unbound.
    """
    N = len(bound)
    block_size = N // n_blocks
    deltaGs, p_b, p_ub = [], [], []
    for i in range(n_blocks):
        start, end = i * block_size, (i + 1) * block_size if i < n_blocks - 1 else N
        block = bound[start:end]
        p_bound = block.mean()
        p_unbound = 1 - p_bound
        if p_bound > 0 and p_unbound > 0:
            deltaG = -kT * np.log(p_bound / p_unbound)
            deltaGs.append(deltaG)
            p_b.append(p_bound)
            p_ub.append(p_unbound)
        else:
            deltaGs.append(np.nan)
            p_b.append(np.nan)
            p_ub.append(np.nan)
    deltaGs, p_b, p_ub = map(np.array, (deltaGs, p_b, p_ub))
    valid = ~np.isnan(deltaGs)
    if valid.sum() > 0:
        stats = (
            deltaGs[valid].mean(),
            deltaGs[valid].std(ddof=1),
            p_b[valid].mean(),
            p_b[valid].std(ddof=1),
            p_ub[valid].mean(),
            p_ub[valid].std(ddof=1),
        )
    else:
        stats = (None,) * 6
    return stats


def analyze_binding_trajectory(
    traj_file,
    top_pdb_file,
    binding_site,
    ligand_resname,
    dist_cutoff=1.5,
    n_blocks=5,
    centroid_cutoff=0.25,
    output_bound="bound_unbound_states.npy",
):
    """
    Analyze the trajectory to determine bound/unbound states and estimate ΔG.
    - Loads trajectory and topology
    - Applies atom filtering and centroid calculation
    - Saves bound/unbound state array
    - Returns (stats, bound_array, centroid_distances)
    """
    pdb = PDBFile(top_pdb_file)
    openmm_top = pdb.topology
    md_top = md.load_topology(top_pdb_file)
    traj = md.load(traj_file, top=md_top)
    site_indices = get_binding_site_indices(openmm_top, binding_site)
    ligand_indices = get_ligand_indices(openmm_top, ligand_resname)
    ligand_indices_filt = filter_atoms_near_centroid(
        traj.xyz[0], ligand_indices, centroid_cutoff
    )
    site_indices_filt = filter_atoms_near_centroid(
        traj.xyz[0], site_indices, centroid_cutoff
    )
    overlap = set(ligand_indices_filt) & set(site_indices_filt)
    ligand_indices_filt = [i for i in ligand_indices_filt if i not in overlap]
    site_indices_filt = [i for i in site_indices_filt if i not in overlap]

    if not ligand_indices_filt or not site_indices_filt:
        print("Error: CV atom selection failed in analysis.")
        sys.exit(1)

    ligand_coords = traj.xyz[:, ligand_indices_filt, :].mean(axis=1)
    site_coords = traj.xyz[:, site_indices_filt, :].mean(axis=1)
    centroid_dist = np.linalg.norm(ligand_coords - site_coords, axis=1)
    bound = centroid_dist < dist_cutoff
    np.save(output_bound, bound)  # Save bound/unbound state for further analysis

    stats = block_average_deltaG(bound, kT=2.479, n_blocks=n_blocks)
    return stats, bound, centroid_dist


def main():
    """
    Main entry point for the pipeline.
    - Parses arguments
    - Loads system and binding site
    - Selects atoms for CV
    - Runs MD and metadynamics
    - Analyzes trajectory and saves results
    """
    parser = argparse.ArgumentParser(
        description="Run ligand-DNA binding/unbinding metadynamics and analyze states."
    )
    parser.add_argument("--system", required=True, help="OpenMM system.xml")
    parser.add_argument("--pdb", required=True, help="Input structure PDB")
    parser.add_argument("--site", required=True, help="Binding site JSON file")
    parser.add_argument(
        "--ligand", required=True, help="Ligand residue name (e.g., X8V)"
    )
    parser.add_argument("--steps", type=int, default=500000, help=" 500k ~ 1 ns")
    parser.add_argument(
        "--interval", type=int, default=1000, help="Trajectory/report interval"
    )
    parser.add_argument("--traj", default="traj.dcd", help="Output trajectory file")
    parser.add_argument(
        "--dist_cutoff", type=float, default=1.0, help="Centroid distance cutoff in nm"
    )
    parser.add_argument(
        "--blocks", type=int, default=5, help="Number of blocks for error estimation"
    )
    parser.add_argument(
        "--centroid_cutoff",
        type=float,
        default=0.25,
        help="Max distance from centroid for CV group (nm)",
    )
    parser.add_argument(
        "--metad_height",
        type=float,
        default=15,
        help="Metadynamics Gaussian height (kJ/mol)",
    )
    parser.add_argument(
        "--metad_sigma_dist",
        type=float,
        default=0.3,
        help="Metadynamics Gaussian width (nm)",
    )
    parser.add_argument(
        "--metad_biasfactor", type=float, default=10.0, help="Metadynamics bias factor"
    )
    parser.add_argument(
        "--metad_stride",
        type=int,
        default=1000,
        help="Metadynamics hill deposition stride",
    )
    parser.add_argument(
        "--wall_buffer", type=float, default=2.5, help="Buffer (nm) for wall restraint"
    )
    parser.add_argument(
        "--k_wall",
        type=float,
        default=20.0,
        help="Force constant for wall restraint (kJ/mol/nm^2)",
    )
    parser.add_argument(
        "--nvt_steps", type=int, default=10000, help="NVT equilibration steps"
    )
    parser.add_argument(
        "--npt_steps", type=int, default=20000, help="NPT equilibration steps"
    )
    parser.add_argument(
        "--eq_steps", type=int, default=10000, help="Final NVT equilibration steps"
    )
    parser.add_argument(
        "--temperature", type=float, default=300, help="Simulation temperature (K)"
    )
    parser.add_argument(
        "--pressure", type=float, default=1.0, help="Pressure for NPT (bar)"
    )
    parser.add_argument(
        "--output_prefix", default="output", help="Prefix for output files"
    )
    args = parser.parse_args()

    t0 = time.time()
    pdb = PDBFile(args.pdb)
    topology = pdb.topology
    binding_site = load_binding_site(args.site)
    site_indices = get_binding_site_indices(topology, binding_site)
    ligand_indices = get_ligand_indices(topology, args.ligand)

    # Filter CV atoms to those near the centroid (removes flexible tails, etc.)
    ligand_indices_filt = filter_atoms_near_centroid(
        pdb.positions, ligand_indices, args.centroid_cutoff
    )
    site_indices_filt = filter_atoms_near_centroid(
        pdb.positions, site_indices, args.centroid_cutoff
    )
    # Remove any overlap between site and ligand indices
    overlap = set(ligand_indices_filt) & set(site_indices_filt)
    ligand_indices_filt = [i for i in ligand_indices_filt if i not in overlap]
    site_indices_filt = [i for i in site_indices_filt if i not in overlap]

    if not ligand_indices_filt or not site_indices_filt:
        print("Error: CV atom selection failed.")
        sys.exit(1)

    traj_file, _, equilibrated_pdb = run_md_with_metadynamics(
        args.system,
        args.pdb,
        ligand_indices_filt,
        site_indices_filt,
        nsteps=args.steps,
        report_interval=args.interval,
        output_traj=args.traj,
        metad_height=args.metad_height,
        metad_sigma_dist=args.metad_sigma_dist,
        metad_biasfactor=args.metad_biasfactor,
        metad_stride=args.metad_stride,
        temperature=args.temperature,
        pressure=args.pressure,
        wall_buffer=args.wall_buffer,
        k_wall=args.k_wall,
        nvt_steps=args.nvt_steps,
        npt_steps=args.npt_steps,
        eq_steps=args.eq_steps,
        output_prefix=args.output_prefix,
    )

    stats, _, _ = analyze_binding_trajectory(
        traj_file,
        equilibrated_pdb,
        binding_site,
        args.ligand,
        dist_cutoff=args.dist_cutoff,
        n_blocks=args.blocks,
        centroid_cutoff=args.centroid_cutoff,
        output_bound=f"{args.output_prefix}_bound_unbound.npy",
    )

    t1 = time.time()
    elapsed_sec = t1 - t0
    deltaG, deltaG_err, p_bound, p_bound_err, p_unbound, p_unbound_err = stats
    with open("binding_energy_result.json", "w") as f:
        json.dump(
            {
                "deltaG_kJ_per_mol": deltaG,
                "deltaG_error_kJ_per_mol": deltaG_err,
                "P_bound": p_bound,
                "P_bound_error": p_bound_err,
                "P_unbound": p_unbound,
                "P_unbound_error": p_unbound_err,
                "wall_time_sec": elapsed_sec,
            },
            f,
            indent=2,
        )
    print(f"Bound/unbound state saved to {args.output_prefix}_bound_unbound.npy")
    print(
        f"Centroid distances saved for analysis. ΔG = {deltaG:.2f} ± {deltaG_err:.2f} kJ/mol"
    )


if __name__ == "__main__":
    main()
