# run_free_energy.py

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
)
from openmm.unit import kelvin, bar
import mdtraj as md
from openmm.app.metadynamics import Metadynamics, BiasVariable

SOLVENT_IONS = {"HOH", "WAT", "NA", "CL", "K", "MG", "CA", "SO4", "PO4"}


def strip_unit(x):
    try:
        return np.array(x.value_in_unit(unit.nanometer))
    except AttributeError:
        return np.array(x)


def load_binding_site(json_file):
    with open(json_file, "r") as f:
        return json.load(f)


def is_heavy_atom(atom):
    return atom.element.symbol != "H"


def get_binding_site_indices(topology, binding_site):
    if isinstance(topology, md.Topology):
        md_top = topology
    else:
        md_top = md.Topology.from_openmm(topology)
    table, _ = md_top.to_dataframe()
    pdbid_to_chain_index = {
        getattr(chain, "chain_id", None): i
        for i, chain in enumerate(md_top.chains)
        if getattr(chain, "chain_id", None) is not None
    }
    indices = []
    for res in binding_site:
        if res["resname"].strip().upper() in SOLVENT_IONS:
            continue
        chain_id = res["chain"]
        chain_index = pdbid_to_chain_index.get(chain_id)
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
    if hasattr(topology, "atoms") and not callable(topology.atoms):
        atoms = topology.atoms
    else:
        atoms = topology.atoms()
    return [
        atom.index
        for atom in atoms
        if atom.residue.name == ligand_resname and atom.element.symbol != "H"
    ]


def get_centroid(positions, indices):
    coords = np.array([strip_unit(positions[i]) for i in indices])
    return coords.mean(axis=0)


def filter_atoms_near_centroid(positions, indices, cutoff_nm):
    if len(indices) == 0:
        return []
    centroid = get_centroid(positions, indices)
    coords = np.array([strip_unit(positions[i]) for i in indices])
    dists = np.linalg.norm(coords - centroid, axis=1)
    return [idx for idx, d in zip(indices, dists) if d <= cutoff_nm]


def create_centroid_distance_cv(ligand_indices, site_indices):
    force = CustomCentroidBondForce(2, "distance(g1,g2)")
    force.addGroup(ligand_indices)
    force.addGroup(site_indices)
    force.addBond([0, 1], [])
    return force


def create_ligand_wall_restraint(
    ligand_indices, box_vectors, wall_buffer=2.0, k_wall=100.0
):
    box_lengths = [box_vectors[i][i]._value for i in range(3)]
    lower = [wall_buffer, wall_buffer, wall_buffer]
    upper = [box_lengths[i] - wall_buffer for i in range(3)]
    print(
        f"Wall restraint: lower bounds = {lower}, upper bounds = {upper}, box = {box_lengths}"
    )
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
    force.addBond(
        [0],
        [
            lower[0],
            upper[0],
            lower[1],
            upper[1],
            lower[2],
            upper[2],
            k_wall,
        ],
    )
    return force


def create_centroid_harmonic_restraint(
    ligand_indices, site_indices, k_rest=1000.0, r0=0.8
):
    # Harmonic restraint on the centroid distance, only acts when ligand is far from site
    expr = "0.5 * k_rest * (distance(g1,g2) - r0)^2"
    force = CustomCentroidBondForce(2, expr)
    force.addGroup(ligand_indices)
    force.addGroup(site_indices)
    force.addPerBondParameter("k_rest")
    force.addPerBondParameter("r0")
    force.addBond([0, 1], [k_rest, r0])
    return force


def run_md_metadynamics(
    system_xml,
    pdb_file,
    ligand_indices,
    site_indices,
    nsteps=500000,
    report_interval=500,
    output_traj="traj.dcd",
    metad_height=50.0,
    metad_sigma_dist=0.2,
    metad_biasfactor=5,
    metad_stride=50,
    temperature=300,
    pressure=1.0,
    wall_buffer=2.0,
    k_wall=100.0,
    nvt_steps=5000,
    npt_steps=10000,
    eq_steps=5000,
    output_prefix="output",
    k_rest=1000.0,
    r0=0.8,
):
    # --- Load system and structure ---
    system = XmlSerializer.deserialize(open(system_xml).read())
    pdb = PDBFile(pdb_file)
    box_vectors = pdb.topology.getPeriodicBoxVectors()
    print(f"Initial box vectors (nm): {box_vectors}")

    # --- Add wall restraint ---
    wall_force = create_ligand_wall_restraint(
        ligand_indices,
        box_vectors,
        wall_buffer=wall_buffer,
        k_wall=k_wall,
    )
    system.addForce(wall_force)

    # --- Add harmonic centroid restraint ---
    harmonic_force = create_centroid_harmonic_restraint(
        ligand_indices, site_indices, k_rest=k_rest, r0=r0
    )
    system.addForce(harmonic_force)

    # --- Step 1: Minimization ---
    integrator_nvt = LangevinIntegrator(
        temperature * kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = Simulation(pdb.topology, system, integrator_nvt)
    simulation.context.setPositions(pdb.positions)
    print("Minimizing energy for 5,000 steps...")
    simulation.minimizeEnergy(maxIterations=5000)

    # --- Step 2: NVT Equilibration (short) ---
    simulation.context.setVelocitiesToTemperature(temperature * kelvin)
    print(f"NVT equilibration ({nvt_steps} steps)...")
    simulation.step(nvt_steps)

    # --- Step 3: NPT Equilibration (allow box to relax) ---
    barostat = MonteCarloBarostat(pressure * bar, temperature * kelvin, 25)
    system.addForce(barostat)
    integrator_npt = LangevinIntegrator(
        temperature * kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation_npt = Simulation(pdb.topology, system, integrator_npt)
    simulation_npt.context.setPositions(
        simulation.context.getState(getPositions=True).getPositions()
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
    forces = [system.getForce(i) for i in range(system.getNumForces())]
    for i, force in enumerate(forces):
        if isinstance(force, MonteCarloBarostat):
            system.removeForce(i)
            break
    # Get new box vectors and positions
    state = simulation_npt.context.getState(
        getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True
    )
    positions = state.getPositions()
    velocities = state.getVelocities()
    box_vectors = state.getPeriodicBoxVectors()
    print(f"Box after NPT: {box_vectors}")

    # --- Step 4: NVT Equilibration (short, to stabilize after pressure) ---
    integrator_nvt2 = LangevinIntegrator(
        temperature * kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation_nvt2 = Simulation(pdb.topology, system, integrator_nvt2)
    simulation_nvt2.context.setPositions(positions)
    simulation_nvt2.context.setVelocities(velocities)
    a, b, c = box_vectors
    simulation_nvt2.context.setPeriodicBoxVectors(a, b, c)
    print(f"Final NVT equilibration ({eq_steps} steps)...")
    simulation_nvt2.step(eq_steps)

    # === Save equilibrated structure and state ===
    eq_state = simulation_nvt2.context.getState(
        getPositions=True, enforcePeriodicBox=True
    )
    eq_positions = eq_state.getPositions()
    eq_box = eq_state.getPeriodicBoxVectors()
    equilibrated_pdb = f"{output_prefix}_equilibrated.pdb"

    # Set box vectors in topology before writing
    pdb.topology.setPeriodicBoxVectors(eq_box)
    with open(equilibrated_pdb, "w") as f:
        PDBFile.writeFile(pdb.topology, eq_positions, f, keepIds=True)
    print(f"Saved equilibrated structure to {equilibrated_pdb}")

    # --- Step 5: Production with Metadynamics (NVT) ---
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
    integrator_prod = LangevinIntegrator(
        temperature * kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation_prod = Simulation(pdb.topology, system, integrator_prod)
    simulation_prod.context.setPositions(
        simulation_nvt2.context.getState(getPositions=True).getPositions()
    )
    simulation_prod.context.setVelocities(
        simulation_nvt2.context.getState(getVelocities=True).getVelocities()
    )
    a, b, c = box_vectors
    simulation_prod.context.setPeriodicBoxVectors(a, b, c)
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

    return output_traj, pdb.topology


def block_average_deltaG(bound, kT=2.479, n_blocks=5):
    N = len(bound)
    block_size = N // n_blocks
    deltaGs = []
    p_b = []
    p_ub = []
    for i in range(n_blocks):
        start = i * block_size
        end = (i + 1) * block_size if i < n_blocks - 1 else N
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
        mean_deltaG = std_deltaG = mean_p_bound = std_p_bound = mean_p_unbound = (
            std_p_unbound
        ) = None
    return (
        mean_deltaG,
        std_deltaG,
        mean_p_bound,
        std_p_bound,
        mean_p_unbound,
        std_p_unbound,
    )


def compute_binding_free_energy(
    traj_file,
    top_pdb_file,
    binding_site,
    ligand_resname,
    dist_cutoff=1.5,
    n_blocks=5,
    centroid_cutoff=0.25,
):
    top = md.load_topology(top_pdb_file)
    traj = md.load(traj_file, top=top)
    site_indices = get_binding_site_indices(traj.topology, binding_site)
    ligand_indices = get_ligand_indices(traj.topology, ligand_resname)

    ligand_indices_filt = filter_atoms_near_centroid(
        traj.xyz[0], ligand_indices, centroid_cutoff
    )
    site_indices_filt = filter_atoms_near_centroid(
        traj.xyz[0], site_indices, centroid_cutoff
    )
    overlap = set(ligand_indices_filt) & set(site_indices_filt)
    ligand_indices_filt = [i for i in ligand_indices_filt if i not in overlap]
    site_indices_filt = [i for i in site_indices_filt if i not in overlap]
    if len(ligand_indices_filt) == 0 or len(site_indices_filt) == 0:
        print("Error: CV atom selection failed in analysis.")
        sys.exit(1)

    ligand_coords = traj.xyz[:, ligand_indices_filt, :].mean(axis=1)
    site_coords = traj.xyz[:, site_indices_filt, :].mean(axis=1)
    centroid_dist = np.linalg.norm(ligand_coords - site_coords, axis=1)
    bound = centroid_dist < dist_cutoff
    kT = 2.479
    (
        mean_deltaG,
        std_deltaG,
        mean_p_bound,
        std_p_bound,
        mean_p_unbound,
        std_p_unbound,
    ) = block_average_deltaG(bound, kT=kT, n_blocks=n_blocks)
    return (
        mean_deltaG,
        std_deltaG,
        mean_p_bound,
        std_p_bound,
        mean_p_unbound,
        std_p_unbound,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run MD and estimate binding free energy using centroid distance CV metadynamics."
    )
    parser.add_argument("--system", required=True, help="OpenMM system.xml")
    parser.add_argument("--pdb", required=True, help="Input structure PDB")
    parser.add_argument("--site", required=True, help="Binding site JSON file")
    parser.add_argument(
        "--ligand", required=True, help="Ligand residue name (e.g., DAN)"
    )
    parser.add_argument("--steps", type=int, default=500000, help="Number of MD steps")
    parser.add_argument(
        "--interval", type=int, default=500, help="Trajectory/report interval"
    )
    parser.add_argument("--traj", default="traj.dcd", help="Output trajectory file")
    parser.add_argument(
        "--dist_cutoff", type=float, default=1.5, help="Centroid distance cutoff in nm"
    )
    parser.add_argument(
        "--blocks", type=int, default=5, help="Number of blocks for error estimation"
    )
    parser.add_argument(
        "--centroid_cutoff",
        type=float,
        default=0.25,
        help="Atom must be within this distance (nm) of centroid to be included",
    )
    parser.add_argument(
        "--metad_height",
        type=float,
        default=50.0,
        help="Metadynamics Gaussian height (kJ/mol)",
    )
    parser.add_argument(
        "--metad_sigma_dist",
        type=float,
        default=0.2,
        help="Metadynamics Gaussian width for centroid distance CV (nm)",
    )
    parser.add_argument(
        "--metad_biasfactor", type=float, default=5.0, help="Metadynamics bias factor"
    )
    parser.add_argument(
        "--metad_stride",
        type=int,
        default=50,
        help="Metadynamics hill deposition stride",
    )
    parser.add_argument(
        "--wall_buffer",
        type=float,
        default=2.0,
        help="Buffer (nm) inside box for ligand wall restraint",
    )
    parser.add_argument(
        "--k_wall",
        type=float,
        default=100.0,
        help="Force constant (kJ/mol/nm^2) for ligand wall restraint",
    )
    parser.add_argument(
        "--nvt_steps", type=int, default=5000, help="NVT equilibration steps"
    )
    parser.add_argument(
        "--npt_steps", type=int, default=10000, help="NPT equilibration steps"
    )
    parser.add_argument(
        "--eq_steps", type=int, default=5000, help="Final NVT equilibration steps"
    )
    parser.add_argument(
        "--temperature", type=float, default=300, help="Simulation temperature (K)"
    )
    parser.add_argument(
        "--pressure", type=float, default=1.0, help="Pressure for NPT (bar)"
    )
    parser.add_argument(
        "--output_prefix",
        default="output",
        help="Prefix for output files (e.g. output_equilibrated.pdb)",
    )
    parser.add_argument(
        "--k_rest",
        type=float,
        default=1000.0,
        help="Force constant (kJ/mol/nm^2) for centroid harmonic restraint",
    )
    parser.add_argument(
        "--r0",
        type=float,
        default=0.8,
        help="Equilibrium centroid distance for harmonic restraint (nm)",
    )
    args = parser.parse_args()

    t0 = time.time()
    pdb = PDBFile(args.pdb)
    topology = pdb.topology
    binding_site = load_binding_site(args.site)
    site_indices = get_binding_site_indices(topology, binding_site)
    ligand_indices = get_ligand_indices(topology, args.ligand)
    print(f"Initial ligand indices: {ligand_indices}")
    print(f"Initial site indices: {site_indices}")

    ligand_indices_filt = filter_atoms_near_centroid(
        pdb.positions, ligand_indices, args.centroid_cutoff
    )
    site_indices_filt = filter_atoms_near_centroid(
        pdb.positions, site_indices, args.centroid_cutoff
    )
    overlap = set(ligand_indices_filt) & set(site_indices_filt)
    ligand_indices_filt = [i for i in ligand_indices_filt if i not in overlap]
    site_indices_filt = [i for i in site_indices_filt if i not in overlap]
    print(f"Ligand CV group: {len(ligand_indices_filt)} atoms")
    print(f"Site CV group: {len(site_indices_filt)} atoms")
    if len(ligand_indices_filt) == 0 or len(site_indices_filt) == 0:
        print("Error: CV atom selection failed.")
        sys.exit(1)
    traj_file, top = run_md_metadynamics(
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
        k_rest=args.k_rest,
        r0=args.r0,
    )
    equilibrated_pdb = f"{args.output_prefix}_equilibrated.pdb"
    deltaG, deltaG_err, p_bound, p_bound_err, p_unbound, p_unbound_err = (
        compute_binding_free_energy(
            traj_file,
            equilibrated_pdb,
            binding_site,
            args.ligand,
            dist_cutoff=args.dist_cutoff,
            n_blocks=args.blocks,
            centroid_cutoff=args.centroid_cutoff,
        )
    )
    t1 = time.time()
    elapsed_sec = t1 - t0
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
