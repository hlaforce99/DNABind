# run_fep.py

import argparse
import numpy as np
import json
import time
from openmm.app import PDBFile, Simulation
from openmm import XmlSerializer, unit
import openmm as mm
import sys
import logging

def setup_logging():
    logging.basicConfig(filename='fep.log', level=logging.INFO,
                        format='%(asctime)s %(levelname)s:%(message)s')
    logging.info("Started run_fep.py")

def get_ligand_atom_indices(topology, ligand_resname):
    """Return atom indices for the ligand."""
    return [atom.index for atom in topology.atoms() if atom.residue.name == ligand_resname]

def scale_ligand_charges(system, ligand_indices, orig_charges, lam):
    """Scale the charges of ligand atoms by lambda."""
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            for idx, atom_idx in enumerate(ligand_indices):
                q, sig, eps = force.getParticleParameters(atom_idx)
                force.setParticleParameters(atom_idx, orig_charges[idx]*lam, sig, eps)
            force.updateParametersInContext = getattr(force, 'updateParametersInContext', lambda context: None)
            return force  # Only one NonbondedForce
    raise RuntimeError("No NonbondedForce found in system.")

def run_fep(system_xml, integrator_xml, pdb_file, ligand_resname, nsteps, n_windows, temperature, output_prefix):
    setup_logging()
    logging.info(f"Parameters: system={system_xml}, integrator={integrator_xml}, pdb={pdb_file}, ligand={ligand_resname}, nsteps={nsteps}, n_windows={n_windows}, T={temperature}")

    # Load system, integrator, and structure
    system = XmlSerializer.deserialize(open(system_xml).read())
    integrator = XmlSerializer.deserialize(open(integrator_xml).read())
    pdb = PDBFile(pdb_file)
    topology = pdb.topology

    # Find ligand atoms
    ligand_atoms = get_ligand_atom_indices(topology, ligand_resname)
    if not ligand_atoms:
        logging.error(f"No ligand atoms found with residue name {ligand_resname}.")
        sys.exit(f"No ligand atoms found with residue name {ligand_resname}.")

    # Get original ligand charges
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            orig_charges = [force.getParticleParameters(idx)[0] for idx in ligand_atoms]
            break
    else:
        logging.error("No NonbondedForce found in system.")
        sys.exit("No NonbondedForce found in system.")

    # Lambda schedule
    lambdas = np.linspace(1.0, 0.0, n_windows)
    window_energies = []
    window_stderrs = []
    wall_times = []

    # Simulation setup
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)

    for i, lam in enumerate(lambdas):
        logging.info(f"Lambda window {i+1}/{n_windows}: lambda={lam:.3f}")
        print(f"Lambda window {i+1}/{n_windows}: lambda={lam:.3f}")

        # Scale ligand charges
        scale_ligand_charges(system, ligand_atoms, orig_charges, lam)
        try:
            system.getForces()[0].updateParametersInContext(simulation.context)
        except Exception:
            pass  # For OpenMM < 8, this may not exist

        # Equilibration (short)
        simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
        simulation.step(int(nsteps * 0.1))

        # Production
        energies = []
        t0 = time.time()
        for step in range(nsteps):
            simulation.step(1)
            if step % max(1, nsteps // 10) == 0:
                print(f"  Step {step}/{nsteps}...", end='\r')
            state = simulation.context.getState(getEnergy=True)
            energies.append(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        t1 = time.time()
        wall_times.append(t1 - t0)
        mean_E = np.mean(energies)
        std_E = np.std(energies) / np.sqrt(len(energies))
        window_energies.append(mean_E)
        window_stderrs.append(std_E)
        logging.info(f"Window {i}: mean_E={mean_E:.2f} kJ/mol, stderr={std_E:.2f} kJ/mol, wall_time={t1-t0:.1f}s")

    # TI integration (trapezoidal rule)
    dG = 0.0
    for i in range(1, n_windows):
        delta_lambda = lambdas[i] - lambdas[i-1]
        avg_dE = 0.5 * (window_energies[i] + window_energies[i-1])
        dG += avg_dE * delta_lambda
    dG_kcal = dG / 4.184  # Convert kJ/mol to kcal/mol

    # Output results
    result = {
        "lambdas": list(map(float, lambdas)),
        "mean_energies_kJ_per_mol": list(map(float, window_energies)),
        "stderr_kJ_per_mol": list(map(float, window_stderrs)),
        "wall_times_sec": list(map(float, wall_times)),
        "dG_TI_kJ_per_mol": float(dG),
        "dG_TI_kcal_per_mol": float(dG_kcal),
        "nsteps_per_window": nsteps,
        "n_windows": n_windows,
        "ligand_resname": ligand_resname,
        "method": "Coulomb decoupling, simple TI"
    }
    with open(f"{output_prefix}_fep_result.json", "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nFEP/TI complete. ΔG = {dG:.2f} kJ/mol ({dG_kcal:.2f} kcal/mol)")
    print(f"Results saved to {output_prefix}_fep_result.json")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Alchemical FEP/TI: ligand Coulomb decoupling in binding site")
    parser.add_argument("--system", required=True, help="System XML file")
    parser.add_argument("--integrator", required=True, help="Integrator XML file")
    parser.add_argument("--pdb", required=True, help="Structure PDB file")
    parser.add_argument("--ligand", required=True, help="Ligand residue name (e.g. X8V)")
    parser.add_argument("--nsteps", type=int, default=5000, help="MD steps per lambda window (default: 5000)")
    parser.add_argument("--windows", type=int, default=5, help="Number of lambda windows (default: 5)")
    parser.add_argument("--temperature", type=float, default=300, help="Temperature in K (default: 300)")
    parser.add_argument("--output_prefix", type=str, default="fep", help="Prefix for output files")
    args = parser.parse_args()

    run_fep(
        args.system,
        args.integrator,
        args.pdb,
        args.ligand,
        args.nsteps,
        args.windows,
        args.temperature,
        args.output_prefix
    )
