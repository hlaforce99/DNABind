# run_free_energy.py

import argparse
import json
import numpy as np
from openmm.app import PDBFile, Simulation
from openmm import XmlSerializer
from openmm.unit import kelvin
from sys import stdout
import deeptime
import mdtraj as md
import time

def load_binding_site(json_file):
    with open(json_file, 'r') as f:
        return json.load(f)

def get_binding_site_indices(topology, binding_site):
    indices = []
    for atom in topology.atoms():
        for res in binding_site:
            if (atom.residue.chain.id == res['chain'] and
                int(atom.residue.id) == int(res['resid']) and
                atom.residue.name == res['resname']):
                indices.append(atom.index)
    return indices

def get_ligand_indices(topology, ligand_resname):
    return [atom.index for atom in topology.atoms() if atom.residue.name == ligand_resname]

def run_md(system_xml, integrator_xml, pdb_file, nsteps=500000, report_interval=1000, output_traj='traj.pdb'):
    system = XmlSerializer.deserialize(open(system_xml).read())
    integrator = XmlSerializer.deserialize(open(integrator_xml).read())
    pdb = PDBFile(pdb_file)
    from openmm.app import Simulation, PDBReporter, StateDataReporter
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    print("Minimizing...")
    simulation.minimizeEnergy()
    print("Equilibrating...")
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    simulation.step(5000)
    print(f"Running production MD for {nsteps*0.002/1000:.1f} ns...")
    simulation.reporters.append(PDBReporter(output_traj, report_interval))
    simulation.reporters.append(StateDataReporter(stdout, report_interval, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=nsteps, separator='\t'))
    simulation.step(nsteps)
    print("MD complete.")
    return output_traj, pdb.topology

def block_average_deltaG(bound, kT=2.479, n_blocks=5):
    N = len(bound)
    block_size = N // n_blocks
    deltaGs = []
    p_b = []
    p_ub = []
    for i in range(n_blocks):
        start = i * block_size
        end = (i+1) * block_size if i < n_blocks-1 else N
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
    # Remove NaNs for stats
    deltaGs = np.array(deltaGs)
    p_b = np.array(p_b)
    p_ub = np.array(p_ub)
    valid = ~np.isnan(deltaGs)
    if valid.sum() > 0:
        mean_deltaG = deltaGs[valid].mean()
        std_deltaG = deltaGs[valid].std(ddof=1)
        mean_p_bound = p_b[valid].mean()
        std_p_bound = p_b[valid].std(ddof=1)
        mean_p_unbound = p_ub[valid].mean()
        std_p_unbound = p_ub[valid].std(ddof=1)
    else:
        mean_deltaG = std_deltaG = mean_p_bound = std_p_bound = mean_p_unbound = std_p_unbound = None
    return mean_deltaG, std_deltaG, mean_p_bound, std_p_bound, mean_p_unbound, std_p_unbound

def compute_binding_free_energy(traj_file, topology, binding_site, ligand_resname, cutoff=0.5, n_blocks=5):
    traj = md.load(traj_file, top=topology)
    site_indices = get_binding_site_indices(traj.topology, binding_site)
    ligand_indices = get_ligand_indices(traj.topology, ligand_resname)
    pairs = np.array([[i, j] for i in ligand_indices for j in site_indices])
    distances = md.compute_distances(traj, pairs)
    min_dist = distances.min(axis=1)
    bound = min_dist < cutoff
    kT = 2.479  # kJ/mol at 298K

    # Block averaging
    mean_deltaG, std_deltaG, mean_p_bound, std_p_bound, mean_p_unbound, std_p_unbound = block_average_deltaG(bound, kT=kT, n_blocks=n_blocks)

    if mean_deltaG is not None:
        print(f"Estimated binding free energy: ΔG = {mean_deltaG:.2f} ± {std_deltaG:.2f} kJ/mol (block average, {n_blocks} blocks)")
        print(f"P_bound = {mean_p_bound:.3f} ± {std_p_bound:.3f}")
        print(f"P_unbound = {mean_p_unbound:.3f} ± {std_p_unbound:.3f}")
    else:
        print("Insufficient transitions between bound/unbound states for ΔG estimate.")
    return mean_deltaG, std_deltaG, mean_p_bound, std_p_bound, mean_p_unbound, std_p_unbound

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MD and estimate binding free energy using deeptime.")
    parser.add_argument("--system", required=True, help="OpenMM system.xml")
    parser.add_argument("--integrator", required=True, help="OpenMM integrator.xml")
    parser.add_argument("--pdb", required=True, help="Input structure PDB")
    parser.add_argument("--site", required=True, help="Binding site JSON file")
    parser.add_argument("--ligand", required=True, help="Ligand residue name (e.g., DAN)")
    parser.add_argument("--steps", type=int, default=500000, help="Number of MD steps (default: 500k ~ 1ns)")
    parser.add_argument("--interval", type=int, default=1000, help="Trajectory/report interval")
    parser.add_argument("--traj", default="traj.pdb", help="Output trajectory file")
    parser.add_argument("--cutoff", type=float, default=0.5, help="Binding cutoff in nm (default: 0.5 nm = 5 Å)")
    parser.add_argument("--blocks", type=int, default=5, help="Number of blocks for error estimation (default: 5)")
    args = parser.parse_args()

    t0 = time.time()

    # Run MD
    traj_file, top = run_md(args.system, args.integrator, args.pdb, nsteps=args.steps, report_interval=args.interval, output_traj=args.traj)

    # Load binding site
    binding_site = load_binding_site(args.site)

    # Estimate binding free energy with error bars
    deltaG, deltaG_err, p_bound, p_bound_err, p_unbound, p_unbound_err = compute_binding_free_energy(
        traj_file, top, binding_site, args.ligand, cutoff=args.cutoff, n_blocks=args.blocks
    )

    t1 = time.time()
    elapsed_sec = t1 - t0
    print(f"Total wall time: {elapsed_sec/60:.1f} min")

    # Save results
    with open("binding_energy_result.json", "w") as f:
        json.dump({
            "deltaG_kJ_per_mol": deltaG,
            "deltaG_error_kJ_per_mol": deltaG_err,
            "P_bound": p_bound,
            "P_bound_error": p_bound_err,
            "P_unbound": p_unbound,
            "P_unbound_error": p_unbound_err,
            "wall_time_sec": elapsed_sec
        }, f, indent=2)
    print("Results saved to binding_energy_result.json")
