#!/bin/bash

# Usage: ./wrapper.sh PDB_ID LIGAND_RESNAME [BOX_SIZE_NM] [MD_STEPS] [BLOCKS]
# Example: ./wrapper.sh 7KWK X8V 1.0 1000000 5

set -e

if [ $# -lt 2 ]; then
  echo "Usage: $0 PDB_ID LIGAND_RESNAME [BOX_SIZE_NM] [MD_STEPS] [BLOCKS]"
  exit 1
fi

PDBID=$1
LIGAND=$2
BOX=${3:-1.0}
STEPS=${4:-1000000}
BLOCKS=${5:-5}

PREFIX="output"

echo "=== [1/3] Building system with build.py ==="
python build.py --pdb "$PDBID" --prefix "$PREFIX" --box "$BOX" --ligand "$LIGAND"

echo "=== [2/3] Defining binding site with define_binding_site.py ==="

# --- Automatic mode (default): residues within 5 Å of ligand ---
python define_binding_site.py --pdb "${PREFIX}_structure.pdb" --ligand "$LIGAND" --cutoff 5.0 --out binding_site.json

# --- Manual mode: uncomment and edit the next line to specify residues manually ---
# python define_binding_site.py --pdb "${PREFIX}_structure.pdb" --residues "A:10" "B:15" --out binding_site.json

echo "=== [3/3] Running MD and free energy calculation with run_free_energy.py ==="
python run_free_energy.py \
  --system "${PREFIX}_system.xml" \
  --pdb "${PREFIX}_structure.pdb" \
  --site binding_site.json \
  --ligand "$LIGAND" \
  --steps "$STEPS" \
  --interval 1000 \
  --traj traj.dcd \
  --blocks "$BLOCKS" \
  --output_prefix "${PREFIX}"

# === (Optional) ML Potential Integration ===
# To use a TorchMD-NET ML potential, set ML_MODEL_PATH and uncomment the following:
# if [ ! -z "$ML_MODEL_PATH" ]; then
#   python build.py --pdb "$PDBID" --prefix "$PREFIX" --box "$BOX" --ligand "$LIGAND" --ml_potential --ml_model_path "$ML_MODEL_PATH"
# fi

echo "=== Pipeline complete! ==="
echo "Key outputs:"
echo "  - ${PREFIX}_structure.pdb"
echo "  - ${PREFIX}_system.xml"
echo "  - ${PREFIX}_integrator.xml"
echo "  - binding_site.json"
echo "  - traj.dcd"
echo "  - binding_energy_result.json (see this for your binding free energy estimate with error bars and timing)"
