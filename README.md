# DNABind: DNA/RNA-Ligand Binding Free Energy Pipeline

This pipeline estimates the binding free energy of a small molecule ligand to a DNA or RNA target using OpenMM molecular dynamics and a population-based estimator (with deeptime and MDTraj).  
It's fully automated and reproducible: just provide a PDB ID and ligand residue name.

---

## **Requirements**

- [conda](https://docs.conda.io/)
- Python 3.9+
- Packages (install with conda):
  ```bash
  conda install -c conda-forge openmm pdbfixer mdtraj ambertools biopython numpy
  ```
- (Optional, for ML potential integration):
    - [TorchMD-NET](https://github.com/torchmd/torchmd-net?locale=en) (for ML-based force fields)
      ```bash
      pip install torchmd-net
      ```
---

## **Pipeline Overview**

1. **System Building**: Downloads and prepares a solvated, parameterized system from a PDB ID.
2. **Binding Site Definition**:
   - **Automatic**: Identifies residues within a specified distance of the ligand.
   - **Manual**: Lets you specify residues by chain and number.
3. **MD Simulation & Free Energy Estimate**: Runs MD and computes a binding free energy using a simple population estimator, with error bars via block averaging.

---

## **Quick Start**

### **1. Make the wrapper executable**
```bash
chmod +x run.sh
```

### **2. Run the pipeline**
```bash
./run.sh PDB_ID LIGAND_RESNAME [BOX_SIZE_NM] [MD_STEPS] [BLOCKS]
```
- `PDB_ID`: e.g., `7KWK`
- `LIGAND_RESNAME`: 3-letter ligand code, e.g., `X8V`
- `BOX_SIZE_NM` (optional): Water box padding in nm (default: 1.0)
- `MD_STEPS` (optional): Number of MD steps (default: 500000; 500,000 × 2 fs = 1 ns)
- `BLOCKS` (optional): Number of blocks for error estimation (default: 5)

**Example:**
```bash
./run.sh 7KWK X8V
```
or with custom MD steps and blocks:
```bash
./run.sh 7KWK X8V 1.0 1000000 10
```
---

## **ML Potential Integration (TorchMD-NET)**

You can optionally replace the classical bonded and van der Waals (vdW) force fields terms with a machine learning potential trained with [TorchMD-NET](https://github.com/torchmd/torchmd-net?locale=en).

- **When enabled, the pipeline will:**
    - Remove all classical bonded and vdW terms from the OpenMM system.
    - Add a TorchMD-NET force using your trained model (`model.pt`).
      
- **Requirements:**
    - Install TorchMD-NET:
    ```bash
    pip install torchmd-net
    ```
    - Provide a trained TorchMD-NET model checkpoint (see below).

- **How to enable:**
    - Add `--ml_potential` and `--ml_model_path /path/to/model.pt` when calling `build.py` (or modify `wrapper.sh` to add these arguments).
**Example:**
```bash
python build.py --pdb 7KWK --prefix 7KWK_system --box 1.0 --ml_potential --ml_model_path /path/to/best_model.pt
```
Or, if you add these variables to your `wrapper.sh`, you can run:
```bash
ML_MODEL_PATH=/path/to/best_model.pt ./wrapper.sh 7KWK X8V
```

- **How to train a TorchMD-NET model:**
    - Prepare a dataset of structures with energies/forces (see [TorchMD-NET docs](https://torchmd-net.readthedocs.io/en/latest/?locale=en)).
    - Write a config YAML file specifying the training parameters and data paths.
    - Train with:
      ```bash
      torchmd-train --config config.yaml
      ```
    - Use the resulting `best_model.pt` as your `--ml_model_path`.
---
---

## **Manual Binding Site Selection**

If you want to specify the binding site manually, edit the wrapper script or run the step directly:

**Edit the wrapper script** to comment out the automatic line and add your manual selection, e.g.:
```bash
# python binding_site.py --pdb <PDBID>_system_structure.pdb --ligand <LIGAND> --cutoff 5.0 --out binding_site.json
python binding_site.py --pdb <PDBID>_system_structure.pdb --residues "A:10" "A:12" "B:5" --out binding_site.json
```
Or run this command manually before continuing the pipeline.

---

## **Outputs**

- `<PDBID>_system_structure.pdb`: Solvated, fixed structure
- `<PDBID>_system_system.xml`: OpenMM system file
- `<PDBID>_system_integrator.xml`: OpenMM integrator file
- `binding_site.json`: Binding site residue list (automatic or manual)
- `traj.pdb`: MD trajectory (PDB format)
- `binding_energy_result.json`: Binding free energy estimate, with error bars and wall time

---

## **File Structure**

```
├── src/
│   ├── analysis.py
│   ├── binding_site.py
│   ├── build.py
│   ├── config.py
│   ├── main.py
|   ├── simulation.py
|   ├── utils.py
├── README.md
├── environment.yml (optional)
├── run.sh

```

---

## **How It Works**

- **System Preparation:** Uses AMBER force fields with OpenMM and PDBFixer.
- **Binding Site:**
  - **Automatic:** Residues within 5 Å of ligand atoms.
  - **Manual:** User-specified residues.
- **Simulation:** Minimization, equilibration, and production MD (default 1 ns).
- **Free Energy:** Classifies each frame as "bound" or "unbound" and estimates  
  $\[
  \Delta G = -kT \ln\left(\frac{P_\text{bound}}{P_\text{unbound}}\right)
  \]$
  - **Error bars** are computed by block averaging over the trajectory (number of blocks set by `--blocks`).

---

## **Reproducibility**

- All steps are scripted and deterministic.
- Input files and parameters are recorded in the output files.
- To reproduce a result, rerun the pipeline with the same arguments.

---

## **References**

- OpenMM: [http://openmm.org/](http://openmm.org/)
- PDBFixer: [https://github.com/openmm/pdbfixer](https://github.com/openmm/pdbfixer)
- MDTraj: [http://mdtraj.org/](http://mdtraj.org/)
- TorchMD-NET: [https://github.com/torchmd/torchmd-net](https://github.com/torchmd/torchmd-net)

---

## **Contact**

For questions or issues, contact Hunter La Force at hunterlaforce61@gmail.com.
