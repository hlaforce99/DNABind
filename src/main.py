#!/usr/bin/env python3
"""
DNABind 2.0

A pipeline for calculating DNA-ligand binding free energies using metadynamics.

This script integrates all components of the DNABind pipeline:
1. System preparation and building
2. Binding site definition
3. Metadynamics simulation
4. Trajectory analysis and free energy calculation
"""

import os
import sys
import time
import json

# Import DNABind modules
from config import parse_command_line
from utils import (
    check_dependencies,
    cleanup_files,
    create_directory,
    to_serializable,
)
from build import build_system
from binding_site import define_binding_site
from simulation import setup_and_run_metadynamics
from analysis import analyze_trajectory, plot_pmf


def run_dnabind_pipeline(config):
    """
    Run the complete DNABind pipeline with the given configuration.

    Args:
        config: DNABindConfig object

    Returns:
        dict: Results including binding free energy
    """
    start_time = time.time()

    # Set output paths and pdb_id
    output_dir = config.output_dir
    output_prefix = os.path.join(output_dir, config.output_prefix)
    pdb_id = config.pdb_id

    # Create output directory
    create_directory(output_dir)

    # Save configuration
    config.save(os.path.join(output_dir, "dnabind_config.json"))

    # Check dependencies
    if not check_dependencies():
        sys.exit(1)

    print("\n=== System Preparation ===")
    # Build system
    structure_pdb, system_xml = build_system(
        config.pdb_id,
        output_prefix=output_prefix,
        box_padding=config.box_padding,
        ml_potential=config.ml_potential,
        ml_model_path=config.ml_model_path,
        ligand_resname=config.ligand_resname,
        ligand_charge=config.ligand_charge,
        pdb_is_file=config.pdb_is_file,
    )

    print("\n=== Binding Site Definition ===")
    # Define binding site
    binding_site_file = os.path.join(output_dir, "binding_site.json")
    define_binding_site(
        structure_pdb,
        ligand_resname=config.ligand_resname,
        cutoff=config.binding_site_cutoff,
        output_file=binding_site_file,
    )

    print("\n=== Metadynamics Simulation ===")
    # Run metadynamics
    output_traj = os.path.join(output_dir, "traj.dcd")
    traj_file, bias_file, equilibrated_pdb = setup_and_run_metadynamics(
        system_xml,
        structure_pdb,
        binding_site_file,
        config.ligand_resname,
        output_prefix=output_prefix,
        output_traj=output_traj,
        nsteps=config.md_steps,
        report_interval=config.interval,
        metad_height=config.metad_height,
        metad_sigma=config.metad_sigma,
        metad_biasfactor=config.metad_biasfactor,
        metad_stride=config.metad_stride,
        temperature=config.temperature,
        pressure=config.pressure,
        wall_buffer=config.wall_buffer,
        k_wall=config.k_wall,
        nvt_steps=config.nvt_steps,
        npt_steps=config.npt_steps,
        eq_steps=config.eq_steps,
        centroid_cutoff=config.centroid_cutoff,
    )

    print("\n=== Trajectory Analysis ===")
    # Analyze trajectory
    output_bound = os.path.join(output_dir, f"{config.output_prefix}_bound_unbound.npy")
    results = analyze_trajectory(
        traj_file,
        bias_file,
        binding_site_file,
        equilibrated_pdb,
        dist_cutoff=config.dist_cutoff,
        output_bound=output_bound,
        temperature=config.temperature,
        n_blocks=config.blocks,
        output_prefix=output_prefix,
    )

    # Save results
    results_file = os.path.join(
        output_dir, f"{config.output_prefix}_binding_energy.json"
    )
    with open(results_file, "w") as f:
        json.dump(to_serializable(results), f, indent=2)

    # Plot PMF
    plot_pmf(results["pmf_1d"], output_prefix)

    # Clean up temporary files if requested
    if config.cleanup:
        print("\n=== Cleaning up temporary files ===")
        cleanup_files(output_prefix, pdb_id=pdb_id)

    # Print summary
    elapsed_time = time.time() - start_time
    print("\n=== DNABind Complete ===")
    print(f"Total runtime: {elapsed_time/60:.1f} minutes")
    print(
        f"Binding free energy: {results['deltaG_kJ_per_mol']:.2f} ± {results['deltaG_kJ_per_mol_err']:.2f} kJ/mol"
    )
    print(f"Results saved to: {results_file}")

    return results


def main():
    """Main entry point for DNABind pipeline"""
    # Parse command line arguments
    config = parse_command_line()

    # Run pipeline
    try:
        run_dnabind_pipeline(config)
        return 0
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
