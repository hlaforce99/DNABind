#!/usr/bin/env python3
"""
config.py

Configuration management for DNABind pipeline.
"""

import os
import json
import argparse


class DNABindConfig:
    """Configuration manager for DNABind pipeline"""

    def __init__(self, config_file=None, **kwargs):
        # Default parameters
        self.pdb_id = None
        self.ligand_resname = None
        self.box_padding = 2.0
        self.md_steps = 500000
        self.interval = None
        self.metad_stride = None
        self.metad_height = 5.0
        self.metad_sigma = 0.25
        self.metad_biasfactor = 15.0
        self.temperature = 300.0
        self.pressure = 1.0
        self.wall_buffer = 4.0
        self.k_wall = 0.1
        self.binding_site_cutoff = 5.0
        self.dist_cutoff = 1.5
        self.centroid_cutoff = -1
        self.nvt_steps = 10000
        self.npt_steps = 20000
        self.eq_steps = 10000
        self.output_prefix = "output"
        self.cleanup = True
        self.blocks = 5
        self.ml_potential = False
        self.ml_model_path = None
        self.ligand_charge = 0
        self.pdb_is_file = False
        self.output_dir = "."

        # Load from file if provided
        if config_file and os.path.exists(config_file):
            self._load_from_file(config_file)

        # Override with kwargs
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

        # Derived parameters
        if not self.interval:
            self.interval = max(1000, self.md_steps // 1000)
        if not self.metad_stride:
            self.metad_stride = max(100, self.interval // 10)

        # Validate required parameters
        self._validate()

    def _load_from_file(self, config_file):
        """Load configuration from JSON file"""
        with open(config_file, "r") as f:
            config_data = json.load(f)
            for key, value in config_data.items():
                if hasattr(self, key):
                    setattr(self, key, value)

    def _validate(self):
        """Validate configuration parameters"""
        if not self.pdb_id:
            raise ValueError("PDB ID or file path is required")
        if not self.ligand_resname:
            raise ValueError("Ligand residue name is required")

        # Create output directory if needed
        if self.output_dir != "." and not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def save(self, filename="dnabind_config.json"):
        """Save configuration to JSON file"""
        with open(os.path.join(self.output_dir, filename), "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    def to_dict(self):
        """Convert configuration to dictionary"""
        return {k: v for k, v in self.__dict__.items() if not k.startswith("_")}

    def __str__(self):
        """String representation of configuration"""
        return json.dumps(self.to_dict(), indent=2)


def parse_command_line():
    """Parse command line arguments and create configuration"""
    parser = argparse.ArgumentParser(
        description="DNABind: DNA-ligand binding free energy pipeline"
    )
    parser.add_argument("--config", help="Configuration file (JSON)")
    parser.add_argument("--pdb", help="PDB ID or file path")
    parser.add_argument("--ligand", help="Ligand residue name")
    parser.add_argument("--charge", type=int, help="Ligand charge", default=0)
    parser.add_argument("--steps", type=int, help="Number of MD steps")
    parser.add_argument("--output", help="Output prefix")
    parser.add_argument("--outdir", help="Output directory")
    parser.add_argument("--temp", type=float, help="Temperature (K)", default=300.0)
    parser.add_argument("--ml", action="store_true", help="Use ML potential")
    parser.add_argument("--ml_model", help="Path to ML model")
    parser.add_argument(
    "--box_padding", type=float, help="Box padding in nm", default=None
)
    parser.add_argument(
        "--no_cleanup", action="store_true", help="Don't clean up temp files"
    )

    args = parser.parse_args()

    # Create config with file first
    config = DNABindConfig(config_file=args.config)

    # Override with command line arguments
    if args.pdb:
        config.pdb_id = args.pdb
        config.pdb_is_file = os.path.exists(args.pdb)
    if args.ligand:
        config.ligand_resname = args.ligand
    if args.charge is not None:
        config.ligand_charge = args.charge
    if args.steps:
        config.md_steps = args.steps
    if args.box_padding:
        config.box_padding = args.box_padding
    if args.output:
        config.output_prefix = args.output
    if args.outdir:
        config.output_dir = args.outdir
    if args.temp:
        config.temperature = args.temp
    if args.ml:
        config.ml_potential = True
    if args.ml_model:
        config.ml_model_path = args.ml_model
    if args.no_cleanup:
        config.cleanup = False

    return config
