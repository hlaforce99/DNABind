#!/bin/bash

# wrapper.sh
# Usage: ./wrapper.sh PDB_ID LIGAND_RESNAME [BOX_SIZE_NM] [MD_STEPS]

set -e

if [ $# -lt 2 ]; then
  echo "Usage: $0 PDB_ID LIGAND_RESNAME [BOX_SIZE_NM] [MD_STEPS]"
  exit 1
fi

PDBID=$1
LIGAND=$2
BOX=${3:-1.0}
STEPS=${4:-500000}

PREFIX="${PDBID}_system"

echo "=== [1/3] Building system with build.py ==="
python build.py --pdb $PDBID --prefix $PREFIX --box $BOX

echo "=== [2/3] Defining binding site with define_binding_site.py ==="
python define_binding_site.py --pdb ${PREFIX}_structure.pdb --ligand $LIGAND --cutoff 5.0 --out binding_site.json

echo "=== [3/3] Running MD and free energy calculation with run_free_energy.py ==="
python run_free_energy.py \
  --system ${PREFIX}_system.xml \
  --integrator ${PREFIX}_integrator.xml \
  --pdb ${PREFIX}_structure.pdb \
  --site binding_site.json \
  --ligand $LIGAND \
  --steps $STEPS \
  --interval 1000 \
  --traj traj.pdb \
  --cutoff 0.5

echo "=== Pipeline complete! ==="
echo "Check binding_energy_result.json for your binding free energy estimate."

