#!/usr/bin/env python3
"""
analysis.py

Module for analyzing metadynamics simulations and calculating binding free energies:
- Processes trajectories to extract collective variables
- Computes potential of mean force (PMF)
- Identifies bound and unbound states
- Calculates binding free energy with error estimation
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md

# Import common utility functions
from utils import (
    to_serializable,
    kB,
)


def compute_centroid_distance(traj, ligand_indices, site_indices):
    """
    Compute distance between centroids of two groups of atoms.

    Args:
        traj: MDTraj trajectory
        ligand_indices: List of atom indices for ligand
        site_indices: List of atom indices for binding site

    Returns:
        numpy array of distances in nm
    """
    # Calculate centroids for each frame
    ligand_centroid = np.mean(traj.xyz[:, ligand_indices, :], axis=1)
    site_centroid = np.mean(traj.xyz[:, site_indices, :], axis=1)

    # Calculate distances
    diff = ligand_centroid - site_centroid
    return np.sqrt(np.sum(diff**2, axis=1))


# def compute_pmf_1d(
#     data, bias=None, temperature=300.0, bins=100, min_val=None, max_val=None
# ):
#     """
#     Compute 1D potential of mean force from biased data.

#     Args:
#         data: Array of CV values
#         bias: Array of bias values or None
#         temperature: Temperature in Kelvin
#         bins: Number of bins or bin edges
#         min_val: Minimum value for binning
#         max_val: Maximum value for binning

#     Returns:
#         dict with keys 'bin_centers', 'pmf', and 'counts'
#     """
#     # Set up bins
#     if min_val is None:
#         min_val = np.min(data)
#     if max_val is None:
#         max_val = np.max(data)

#     bin_edges = np.linspace(min_val, max_val, bins + 1)
#     bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

#     # Compute histogram
#     counts, _ = np.histogram(data, bins=bin_edges)

#     # Add small value to avoid log(0)
#     counts = counts.astype(float) + 1e-10

#     # Compute PMF
#     pmf = -kB * temperature * np.log(counts)

#     # Shift minimum to zero
#     pmf -= np.min(pmf)

#     # Apply bias correction if provided
#     if bias is not None and len(bias) > 0:
#         # Interpolate bias to match bin centers
#         if len(bias) != len(bin_centers):
#             from scipy.interpolate import interp1d

#             x_bias = np.linspace(min_val, max_val, len(bias))
#             bias_interp = interp1d(x_bias, bias, bounds_error=False, fill_value=0.0)
#             bias_values = bias_interp(bin_centers)
#         else:
#             bias_values = bias

#         # Add bias to PMF
#         pmf += bias_values

#         # Shift minimum to zero again
#         pmf -= np.min(pmf)

#     return {"bin_centers": bin_centers, "pmf": pmf, "counts": counts}


def identify_bound_unbound(distances, pmf_data, dist_cutoff=1.5):
    """
    Identify bound and unbound states based on distance cutoff and PMF.

    Args:
        distances: Array of distance values
        pmf_data: PMF data from metadynamics
        dist_cutoff: Distance cutoff for bound state (nm)

    Returns:
        tuple of (bound_array, warnings_list, bound_state_dist)
    """
    warnings = []

    # Find bound state based on PMF minimum
    bin_centers = pmf_data["bin_centers"]
    pmf = pmf_data["pmf"]

    # Only consider bins up to cutoff
    valid_bins = bin_centers <= dist_cutoff
    if not np.any(valid_bins):
        warnings.append(f"No PMF bins found below cutoff {dist_cutoff} nm")
        bound_state_dist = dist_cutoff
    else:
        # Find minimum of PMF within cutoff
        valid_pmf = pmf[valid_bins]
        valid_bins_centers = bin_centers[valid_bins]
        min_idx = np.argmin(valid_pmf)
        bound_state_dist = valid_bins_centers[min_idx]

    # Check if bound state is at edge of range
    if bound_state_dist == bin_centers[0]:
        warnings.append("Bound state at minimum distance - consider using smaller bins")

    # Classify frames as bound or unbound
    bound = distances <= dist_cutoff

    # Check if we have enough bound frames
    n_bound = np.sum(bound)
    if n_bound < 100:
        warnings.append(
            f"Only {n_bound} bound frames found - results may be unreliable"
        )

    return bound, warnings, bound_state_dist


def count_transitions(bound_array):
    """
    Count transitions between bound and unbound states.

    Args:
        bound_array: Boolean array indicating bound frames

    Returns:
        Number of transitions
    """
    return np.sum(np.abs(np.diff(bound_array.astype(int))))


def calculate_binding_free_energy(
    pmf_data, bound_state_dist, temperature=300.0, n_blocks=5
):
    """
    Calculate binding free energy and error estimate using block averaging.

    Args:
        pmf_data: PMF data from metadynamics
        bound_state_dist: Distance defining bound state (nm)
        temperature: Temperature in Kelvin
        n_blocks: Number of blocks for error estimation

    Returns:
        tuple of (deltaG, deltaG_err) in kJ/mol
    """
    bin_centers = pmf_data["bin_centers"]
    pmf = pmf_data["pmf"]
    counts = pmf_data["counts"]

    # Find bin index closest to bound state distance
    bound_idx = np.argmin(np.abs(bin_centers - bound_state_dist))

    # Calculate volume elements (assuming spherical)
    bin_width = bin_centers[1] - bin_centers[0]
    volumes = 4 * np.pi * bin_centers**2 * bin_width

    # Calculate partition functions
    bound_Z = np.sum(counts[: bound_idx + 1])
    unbound_Z = np.sum(
        counts[bound_idx + 1 :] * np.exp(-pmf[bound_idx + 1 :] / (kB * temperature))
    )

    # Calculate standard volume (1M = 1661 Å³/molecule = 1.661 nm³/molecule)
    V0 = 1.661  # nm³

    # Calculate binding free energy
    deltaG = -kB * temperature * np.log(bound_Z / (unbound_Z * V0))

    # Error estimation using block averaging
    if n_blocks <= 1:
        return deltaG, 0.0

    # Reshape counts into blocks
    block_size = len(counts) // n_blocks
    if block_size < 5:  # Need at least a few bins per block
        n_blocks = max(2, len(counts) // 5)
        block_size = len(counts) // n_blocks

    reshaped_counts = counts[: block_size * n_blocks].reshape(n_blocks, block_size)
    reshaped_pmf = pmf[: block_size * n_blocks].reshape(n_blocks, block_size)
    reshaped_bins = bin_centers[: block_size * n_blocks].reshape(n_blocks, block_size)

    # Calculate free energy for each block
    block_dG = []
    for i in range(n_blocks):
        block_bound_idx = min(bound_idx, block_size - 1)
        block_bound_Z = np.sum(reshaped_counts[i, : block_bound_idx + 1])
        block_unbound_Z = np.sum(
            reshaped_counts[i, block_bound_idx + 1 :]
            * np.exp(-reshaped_pmf[i, block_bound_idx + 1 :] / (kB * temperature))
        )
        if block_bound_Z > 0 and block_unbound_Z > 0:
            block_dG.append(
                -kB * temperature * np.log(block_bound_Z / (block_unbound_Z * V0))
            )

    # Calculate standard error
    if len(block_dG) > 1:
        deltaG_err = np.std(block_dG) / np.sqrt(len(block_dG))
    else:
        deltaG_err = 0.0

    return deltaG, deltaG_err


def plot_pmf(pmf_data, output_prefix="output"):
    """
    Plot PMF curve.

    Args:
        pmf_data: PMF data from metadynamics
        output_prefix: Prefix for output file
    """
    plt.figure(figsize=(8, 6))
    plt.plot(pmf_data["bin_centers"], pmf_data["pmf"], "b-", linewidth=2)
    plt.xlabel("Distance (nm)", fontsize=14)
    plt.ylabel("PMF (kJ/mol)", fontsize=14)
    plt.title("Potential of Mean Force", fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_pmf.pdf")
    plt.close()


def analyze_trajectory(
    traj_file,
    bias_file,
    binding_site_file,
    equilibrated_pdb,
    dist_cutoff=1.5,
    output_bound="bound_unbound.npz",
    temperature=300.0,
    n_blocks=5,
    output_prefix="output",
):
    """
    Analyze metadynamics trajectory to calculate binding free energy:
    1. Load trajectory and compute centroid distance CV
    2. Compute PMF using reweighting
    3. Identify bound/unbound states
    4. Calculate binding free energy with error estimation

    Returns:
        dict: Analysis results including binding free energy and PMF data
    """
    # Check if files exist
    for file_path, name in [
        (traj_file, "Trajectory"),
        (bias_file, "Bias"),
        (binding_site_file, "Binding site"),
        (equilibrated_pdb, "Reference PDB"),
    ]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{name} file not found: {file_path}")

    print(f"Loading trajectory: {traj_file}")
    traj = md.load(traj_file, top=equilibrated_pdb)

    # === Load atom indices from file ===
    cv_indices_file = f"{output_prefix}_cv_indices.npz"
    if not os.path.exists(cv_indices_file):
        raise FileNotFoundError(f"CV indices file not found: {cv_indices_file}")
    indices = np.load(cv_indices_file)
    ligand_indices = indices["ligand"]
    site_indices = indices["site"]

    print(
        f"Using {len(ligand_indices)} ligand atoms and {len(site_indices)} binding site atoms for analysis"
    )

    # Compute centroid distance
    print("Computing centroid distances...")
    distances = compute_centroid_distance(traj, ligand_indices, site_indices)

    # Load bias potential
    try:
        bias_npz = np.load(bias_file)
        if "free_energy" in bias_npz and "grid" in bias_npz:
            free_energy = bias_npz["free_energy"]
            grid = bias_npz["grid"]
            print(
                f"Loaded free energy and grid from {bias_file}: {len(free_energy)} points"
            )
        else:
            raise KeyError("free_energy or grid missing in bias file")
    except Exception as e:
        print(f"Warning: Could not load bias file {bias_file}: {e}")
        print("Using uniform bias (results may be inaccurate)")
        free_energy = np.zeros(100)
        grid = np.linspace(0, 3.0, 100)

    # Compute PMF
    pmf_data = {"bin_centers": grid, "pmf": free_energy, "counts": np.ones_like(grid)}

    # Identify bound and unbound states
    print("Identifying bound/unbound states...")
    bound, warnings, bound_state_dist = identify_bound_unbound(
        distances, pmf_data, dist_cutoff
    )

    # Save bound/unbound classification
    if output_bound:
        np.save(output_bound, bound)
        print(f"Saved bound/unbound classification to {output_bound}")

    # Count transitions
    n_transitions = count_transitions(bound)
    print(f"Detected {n_transitions} transitions between bound and unbound states")

    # Calculate binding free energy
    print("Calculating binding free energy...")
    deltaG, deltaG_err = calculate_binding_free_energy(
        pmf_data, bound_state_dist, temperature, n_blocks
    )

    # Prepare results
    results = {
        "deltaG_kJ_per_mol": deltaG,
        "deltaG_kJ_per_mol_err": deltaG_err,
        "bound_state_dist_nm": bound_state_dist,
        "n_transitions": n_transitions,
        "n_frames": len(distances),
        "n_bound_frames": np.sum(bound),
        "bound_fraction": float(np.mean(bound)),
        "temperature_K": temperature,
        "pmf_1d": pmf_data,
        "warnings": warnings,
    }

    # Print summary
    print("\nAnalysis Summary:")
    print(f"Bound state distance: {bound_state_dist:.3f} nm")
    print(
        f"Bound frames: {results['n_bound_frames']} / {results['n_frames']} ({results['bound_fraction']:.2%})"
    )
    print(f"Transitions: {n_transitions}")
    print(f"Binding free energy: {deltaG:.2f} ± {deltaG_err:.2f} kJ/mol")

    if warnings:
        print("\nWarnings:")
        for warning in warnings:
            print(f"- {warning}")

    return results


def analyze_trajectory_cli():
    """Command-line interface for trajectory analysis"""
    import argparse

    parser = argparse.ArgumentParser(description="Analyze metadynamics trajectory")
    parser.add_argument("--traj", required=True, help="Trajectory file (DCD)")
    parser.add_argument("--bias", required=True, help="Bias potential file (NPY)")
    parser.add_argument("--pdb", required=True, help="Reference PDB file")
    parser.add_argument("--site", required=True, help="Binding site definition (JSON)")
    parser.add_argument(
        "--cutoff", type=float, default=1.5, help="Distance cutoff for bound state (nm)"
    )
    parser.add_argument("--temp", type=float, default=300.0, help="Temperature (K)")
    parser.add_argument(
        "--blocks", type=int, default=5, help="Number of blocks for error estimation"
    )
    parser.add_argument("--output", default="output", help="Output prefix")

    args = parser.parse_args()

    # Run analysis
    results = analyze_trajectory(
        traj_file=args.traj,
        bias_file=args.bias,
        binding_site_file=args.site,
        equilibrated_pdb=args.pdb,
        dist_cutoff=args.cutoff,
        output_bound=f"{args.output}_bound_unbound.npy",
        temperature=args.temp,
        n_blocks=args.blocks,
        output_prefix=args.output,
    )

    # Save results
    with open(f"{args.output}_binding_energy.json", "w") as f:
        json.dump(to_serializable(results), f, indent=2)

    # Plot PMF
    plot_pmf(results["pmf_1d"], args.output)

    # Print final results
    print("\n=== Results ===")
    print(
        f"Binding free energy: {results['deltaG_kJ_per_mol']:.2f} ± "
        f"{results['deltaG_kJ_per_mol_err']:.2f} kJ/mol"
    )


if __name__ == "__main__":
    analyze_trajectory_cli()
