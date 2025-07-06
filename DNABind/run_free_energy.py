# run_free_energy.py

import argparse
import json
import numpy as np
from openmm.app import PDBFile, Simulation
from openmm import XmlSerializer
from openmm.unit import kelvin, picoseconds, nanometers
from sys import stdout
import deeptime
import mdtraj as md

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
    # Load system and integrator
    system = XmlSerializer.deserialize(open(system_xml).read())
    integrator = XmlSerializer.deserialize(open(integrator_xml).read())
    pdb = PDBFile(pdb_file)
    # Set up simulation
    from openmm.app import Simulation, PDBReporter, StateDataReporter
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    # Minimize
    print("Minimizing...")
    simulation.minimizeEnergy()
    # Equilibrate
    print("Equilibrating...")
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    simulation.step(5000)
    # Production
    print(f"Running production MD for {nsteps*0.002/1000:.1f} ns...")
    simulation.reporters.append(PDBReporter(output_traj, report_interval))
    simulation.reporters.append(StateDataReporter(stdout, report_interval, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=nsteps, separator='\t'))
    simulation.step(nsteps)
    print("MD complete.")
    return output_traj, pdb.topology

def compute_binding_free_energy(traj_file, topology, binding_site, ligand_resname, cutoff=0.5):
    # Load trajectory with MDTraj
    traj = md.load(traj_file, top=topology)
    # Get atom indices
    site_indices = get_binding_site_indices(traj.topology, binding_site)
    ligand_indices = get_ligand_indices(traj.topology, ligand_resname)
    # Compute minimum distance between ligand and binding site for each frame
    pairs = np.array([[i, j] for i in ligand_indices for j in site_indices])
    distances = md.compute_distances(traj, pairs)  # shape: (n_frames, n_pairs)
    min_dist = distances.min(axis=1)  # shape: (n_frames,)
    # Bound if min_dist < cutoff (in nm)
    bound = min_dist < cutoff
    # Estimate free energy (ΔG = -kT ln(P_bound/P_unbound))
    kT = 2.479  # kJ/mol at 298K
    p_bound = bound.mean()
    p_unbound = 1 - p_bound
    if p_bound > 0 and p_unbound > 0:
        deltaG = -kT * np.log(p_bound / p_unbound)
        print(f"Estimated binding free energy: ΔG = {deltaG:.2f} kJ/mol")
    else:
        deltaG = None
        print("Insufficient transitions between bound/unbound states for ΔG estimate.")
    return deltaG, p_bound, p_unbound

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
    args = parser.parse_args()

    # Run MD
    traj_file, top = run_md(args.system, args.integrator, args.pdb, nsteps=args.steps, report_interval=args.interval, output_traj=args.traj)

    # Load binding site
    binding_site = load_binding_site(args.site)

    # Estimate binding free energy
    deltaG, p_bound, p_unbound = compute_binding_free_energy(traj_file, top, binding_site, args.ligand, cutoff=args.cutoff)

    # Save results
    with open("binding_energy_result.json", "w") as f:
        json.dump({
            "deltaG_kJ_per_mol": deltaG,
            "P_bound": p_bound,
            "P_unbound": p_unbound
        }, f, indent=2)
    print("Results saved to binding_energy_result.json")

