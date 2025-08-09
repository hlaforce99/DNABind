#!/bin/bash

# DNABind: Wrapper script for DNA-ligand binding free energy calculations
# Usage: ./run.sh [CONFIG_FILE]

set -e

CONFIG_FILE=${1:-"dnabind_config.json"}

# Create default configuration if it doesn't exist
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Creating default configuration file: $CONFIG_FILE"
    cat > "$CONFIG_FILE" <<EOF
{
    "pdb_id": "7KWK",
    "ligand_resname": "X8V",
    "box_padding": 2.0,
    "md_steps": 5000000,
    "metad_height": 15.0,
    "metad_sigma": 0.25,
    "metad_biasfactor": 15.0,
    "temperature": 300.0,
    "pressure": 1.0,
    "binding_site_cutoff": 3.0,
    "centroid_cutoff": 0.5,
    "output_prefix": "output/output",
    "cleanup": true,
    "blocks": 5
}
EOF
    echo "Please edit $CONFIG_FILE and run again."
    exit 1
fi

# Ensure output directory exists
mkdir -p output

# Run the unified pipeline (from src/)
python src/main.py --config "$CONFIG_FILE"

echo "=== Pipeline complete! ==="
echo "Results saved to output/output_binding_energy.json"
