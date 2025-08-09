#!/usr/bin/env python3
"""
simulation.py

Module for running metadynamics simulations with OpenMM:
- System equilibration (NVT and NPT)
- Metadynamics with centroid distance CV
- Trajectory and bias potential output
"""

import sys
import time
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import mdtraj as md
import openmm.app.metadynamics as mtd
from Bio.PDB import PDBParser

from utils import (
    load_binding_site,
    get_binding_site_indices,
    get_ligand_indices,
    filter_atoms_near_centroid,
)


def create_centroid_distance_cv(
    simulation,
    ligand_indices,
    site_indices,
    periodic=True,
    wall_distance=None,
    k_wall=None,
):
    n_atoms = simulation.system.getNumParticles()
    if not ligand_indices or not site_indices:
        raise ValueError("Ligand or binding site has no atoms")
    cv_force = mm.CustomCentroidBondForce(2, "distance(g1, g2)")
    g1 = cv_force.addGroup([i for i in ligand_indices if i < n_atoms])
    g2 = cv_force.addGroup([i for i in site_indices if i < n_atoms])
    cv_force.addBond([g1, g2])
    cv_force.setUsesPeriodicBoundaryConditions(periodic)

    if wall_distance is not None and k_wall is not None:
        energy_expression = (
            f"0.5 * {k_wall} * step(distance(g1, g2) - {wall_distance}) * "
            f"(distance(g1, g2) - {wall_distance})^2"
        )
        wall_force = mm.CustomCentroidBondForce(2, energy_expression)
        g1_wall = wall_force.addGroup([i for i in ligand_indices if i < n_atoms])
        g2_wall = wall_force.addGroup([i for i in site_indices if i < n_atoms])
        wall_force.addBond([g1_wall, g2_wall])
        wall_force.setUsesPeriodicBoundaryConditions(periodic)
        simulation.system.addForce(wall_force)

    # This function is for post-processing, not used during simulation
    def get_cv_value(positions):
        ligand_centroid = np.mean([positions[i] for i in ligand_indices], axis=0)
        site_centroid = np.mean([positions[i] for i in site_indices], axis=0)
        distance = np.linalg.norm(ligand_centroid - site_centroid)
        return distance

    return cv_force, get_cv_value


def add_metadynamics_bias(
    simulation,
    ligand_indices,
    site_indices,
    height=1.0,
    sigma=0.05,
    biasfactor=5.0,
    stride=500,
    temperature=300.0,
    wall_buffer=3.0,
    k_wall=100.0,
    grid_min=0.0,
    grid_max=5.0,
    grid_bins=250,
):
    """
    Add metadynamics bias to simulation.
    """
    simulation.ligand_indices = ligand_indices
    simulation.site_indices = site_indices
    # Compute initial centroid distance for wall placement
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    ligand_centroid = np.mean(
        [positions[i] / unit.nanometer for i in ligand_indices], axis=0
    )
    site_centroid = np.mean(
        [positions[i] / unit.nanometer for i in site_indices], axis=0
    )
    initial_cv = np.linalg.norm(ligand_centroid - site_centroid)
    wall_position = (
        initial_cv + wall_buffer
        if (wall_buffer is not None and k_wall is not None)
        else None
    )
    cv_force, get_cv_value = create_centroid_distance_cv(
        simulation,
        ligand_indices,
        site_indices,
        periodic=True,
        wall_distance=wall_position,
        k_wall=k_wall if wall_position is not None else None,
    )

    meta = mtd.BiasVariable(
        cv_force,
        minValue=grid_min,
        maxValue=grid_max,
        biasWidth=sigma,
        gridWidth=grid_bins,
        periodic=False,
    )
    metadynamics = mtd.Metadynamics(
        system=simulation.system,
        variables=[meta],
        temperature=temperature * unit.kelvin,
        height=height * unit.kilojoule_per_mole,
        biasFactor=biasfactor,
        frequency=stride,
    )

    simulation.context.reinitialize(preserveState=True)
    return simulation, metadynamics, cv_force, get_cv_value, meta


def equilibrate_system(
    system_xml,
    pdb_file,
    output_prefix,
    nvt_steps=10000,
    npt_steps=10000,
    eq_steps=10000,
    temperature=300.0 * unit.kelvin,
    pressure=1.0 * unit.bar,
):
    """
    Equilibrate system with NVT and NPT simulations.
    """
    with open(system_xml, "r") as f:
        system = mm.XmlSerializer.deserialize(f.read())
    pdb = app.PDBFile(pdb_file)
    integrator = mm.LangevinIntegrator(
        temperature, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    print("Minimizing energy...")
    simulation.minimizeEnergy()
    print(f"Running NVT equilibration ({nvt_steps} steps)...")
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(nvt_steps)
    if npt_steps > 0:
        state = simulation.context.getState(getPositions=True)
        positions = state.getPositions()
        system = mm.XmlSerializer.deserialize(open(system_xml).read())
        barostat = mm.MonteCarloBarostat(pressure, temperature)
        system.addForce(barostat)
        integrator = mm.LangevinIntegrator(
            temperature, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
        )
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(positions)
        npt_report_interval = npt_steps // 20
        simulation.reporters.append(
            app.StateDataReporter(
                sys.stdout,
                npt_report_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
            )
        )
        print(f"Running NPT equilibration ({npt_steps} steps)...")
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.step(npt_steps)
    if eq_steps > 0:
        simulation.reporters.clear()
        print(f"Running final equilibration ({eq_steps} steps)...")
        simulation.step(eq_steps)
        state = simulation.context.getState(
            getPositions=True, getVelocities=True, enforcePeriodicBox=True
        )
        positions = state.getPositions()

        # --- Imaging and writing PDB ---
        xyz = (
            np.array(positions.value_in_unit(unit.nanometer))
            .reshape((-1, 3))
            .astype(np.float32)
        )
        top = md.Topology.from_openmm(pdb.topology)
        box_vectors = simulation.context.getState().getPeriodicBoxVectors(asNumpy=True)
        box = np.array(
            [vec.value_in_unit(unit.nanometer) for vec in box_vectors]
        ).astype(np.float32)
        traj = md.Trajectory(xyz=xyz[None, :, :], topology=top)
        traj.unitcell_vectors = box[None, :, :]
        traj.image_molecules(inplace=True)
        equilibrated_pdb = f"{output_prefix}_equilibrated.pdb"
        traj.save_pdb(equilibrated_pdb)
        # --- End imaging section ---

        with open(f"{output_prefix}_equilibrated_system.xml", "w") as f:
            f.write(mm.XmlSerializer.serialize(simulation.system))
        with open(f"{output_prefix}_equilibrated_integrator.xml", "w") as f:
            f.write(mm.XmlSerializer.serialize(simulation.integrator))
        print(f"Equilibration complete. Saved to {equilibrated_pdb}")
        return equilibrated_pdb


def run_metadynamics(
    simulation,
    output_traj="traj.dcd",
    output_bias="bias.npz",
    nsteps=5000000,
    report_interval=5000,
    ligand_indices=None,
    site_indices=None,
    metad_height=1.0,
    metad_sigma=0.05,
    metad_biasfactor=5.0,
    metad_stride=500,
    temperature=300.0,
    wall_buffer=3.0,
    k_wall=100.0,
):
    """
    Run metadynamics simulation with distance-based collective variable.
    """
    simulation.reporters.append(app.DCDReporter(output_traj, report_interval))
    simulation.reporters.append(
        app.StateDataReporter(
            sys.stdout,
            report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=nsteps,
            separator="\t",
        )
    )
    # --- CHANGED: get the BiasVariable object to reconstruct the grid ---
    simulation, metadynamics, _, _, meta = add_metadynamics_bias(
        simulation,
        ligand_indices,
        site_indices,
        height=metad_height,
        sigma=metad_sigma,
        biasfactor=metad_biasfactor,
        stride=metad_stride,
        temperature=temperature,
        wall_buffer=wall_buffer,
        k_wall=k_wall,
    )

    print(f"Running metadynamics for {nsteps} steps...")
    start_time = time.time()
    simulation.step(nsteps)
    end_time = time.time()
    print(f"Simulation complete in {end_time - start_time:.1f} seconds")

    try:
        free_energy = metadynamics.getFreeEnergy()
        # Reconstruct the grid from BiasVariable
        grid = np.linspace(meta.minValue, meta.maxValue, meta.gridWidth)
        np.savez(output_bias, free_energy=free_energy, grid=grid)
        print(f"Saved free energy profile and grid to {output_bias}")
    except AttributeError:
        print("Metadynamics object does not support getFreeEnergy(). No bias saved.")

    return output_traj, output_bias


def setup_and_run_metadynamics(
    system_xml,
    structure_pdb,
    binding_site_file,
    ligand_resname,
    output_prefix="output",
    output_traj="traj.dcd",
    nsteps=5000000,
    report_interval=5000,
    metad_height=1.0,
    metad_sigma=0.05,
    metad_biasfactor=5.0,
    metad_stride=500,
    temperature=300.0,
    pressure=1.0,
    wall_buffer=3.0,
    k_wall=100.0,
    nvt_steps=10000,
    npt_steps=20000,
    eq_steps=10000,
    centroid_cutoff=-1,
):
    # 1. Equilibrate and image system
    equilibrated_pdb = equilibrate_system(
        system_xml,
        structure_pdb,
        output_prefix,
        nvt_steps=nvt_steps,
        npt_steps=npt_steps,
        eq_steps=eq_steps,
        temperature=temperature * unit.kelvin,
        pressure=pressure * unit.bar,
    )

    # 2. Load equilibrated system, positions, and integrator
    with open(f"{output_prefix}_equilibrated_system.xml", "r") as f:
        eq_system = mm.XmlSerializer.deserialize(f.read())
    with open(f"{output_prefix}_equilibrated_integrator.xml", "r") as f:
        eq_integrator = mm.XmlSerializer.deserialize(f.read())
    eq_pdb = app.PDBFile(equilibrated_pdb)
    simulation = app.Simulation(eq_pdb.topology, eq_system, eq_integrator)
    simulation.context.setPositions(eq_pdb.positions)

    # 3. Determine indices
    md_traj = md.load(equilibrated_pdb)
    binding_site = load_binding_site(binding_site_file)
    site_indices = get_binding_site_indices(equilibrated_pdb, binding_site)
    ligand_indices = get_ligand_indices(equilibrated_pdb, ligand_resname)
    if centroid_cutoff > 0:
        site_indices = filter_atoms_near_centroid(
            md_traj.xyz[0], site_indices, centroid_cutoff
        )
        ligand_indices = filter_atoms_near_centroid(
            md_traj.xyz[0], ligand_indices, centroid_cutoff
        )
    overlap = set(ligand_indices) & set(site_indices)
    ligand_indices = [i for i in ligand_indices if i not in overlap]
    site_indices = [i for i in site_indices if i not in overlap]

    # Save indices for analysis
    np.savez(
        f"{output_prefix}_cv_indices.npz", ligand=ligand_indices, site=site_indices
    )
    print(f"Saved CV atom indices to {output_prefix}_cv_indices.npz")

    print(
        f"Using {len(ligand_indices)} ligand atoms and {len(site_indices)} binding site atoms for CV"
    )
    output_bias = f"{output_prefix}_bias.npz"

    # 4. Run metadynamics
    traj_file, bias_file = run_metadynamics(
        simulation,
        output_traj=output_traj,
        output_bias=output_bias,
        nsteps=nsteps,
        report_interval=report_interval,
        ligand_indices=ligand_indices,
        site_indices=site_indices,
        metad_height=metad_height,
        metad_sigma=metad_sigma,
        metad_biasfactor=metad_biasfactor,
        metad_stride=metad_stride,
        temperature=temperature,
        wall_buffer=wall_buffer,
        k_wall=k_wall,
    )
    return traj_file, bias_file, equilibrated_pdb
