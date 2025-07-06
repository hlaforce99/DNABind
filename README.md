# DNA/RNA-Ligand Binding Free Energy Pipeline

This pipeline estimates the binding free energy of a small molecule ligand to a DNA or RNA target using OpenMM molecular dynamics and a population-based estimator (with deeptime and MDTraj).  
It’s fully automated and reproducible: just provide a PDB ID and ligand residue name.

---

## **Requirements**

- [conda](https://docs.conda.io/)
- Python 3.9+
- Packages (install with conda):
  ```bash
  conda install -c conda-forge openmm pdbfixer mdtraj deeptime biopython numpy
  ```

---

## **Pipeline Overview**

1. **System Building**: Downloads and prepares a solvated, parameterized system from a PDB ID.
2. **Binding Site Definition**: Identifies residues near the ligand.
3. **MD Simulation & Free Energy Estimate**: Runs MD and computes a binding free energy using a simple population estimator.

---

## **Quick Start**

### **1. Make the wrapper executable**
```bash
chmod +x wrapper.sh
```

### **2. Run the pipeline**
```bash
./wrapper.sh PDB_ID LIGAND_RESNAME [BOX_SIZE_NM] [MD_STEPS]
```
- `PDB_ID`: e.g., `1D66`
- `LIGAND_RESNAME`: 3-letter ligand code, e.g., `DAN`
- `BOX_SIZE_NM` (optional): Water box padding in nm (default: 1.0)
- `MD_STEPS` (optional): Number of MD steps (default: 500000; 500,000 × 2 fs = 1 ns)

**Example:**
```bash
./wrapper.sh 1D66 DAN
```

---

## **Outputs**

- `<PDBID>_system_structure.pdb`: Solvated, fixed structure
- `<PDBID>_system_system.xml`: OpenMM system file
- `<PDBID>_system_integrator.xml`: OpenMM integrator file
- `binding_site.json`: Binding site residue list
- `traj.pdb`: MD trajectory (PDB format)
- `binding_energy_result.json`: Binding free energy estimate

---

## **File Structure**

```
├── build.py
├── define_binding_site.py
├── run_free_energy.py
├── wrapper.sh
├── README.md
├── environment.yml (optional)
```

---

## **How It Works**

- **System Preparation:** Uses AMBER force fields with OpenMM and PDBFixer.
- **Binding Site:** Residues within 5 Å of ligand atoms.
- **Simulation:** Minimization, equilibration, and production MD (default 1 ns).
- **Free Energy:** Classifies each frame as "bound" or "unbound" and estimates  
ΔG = -kT ln(P_bound / P_unbound)

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
- deeptime: [https://deeptime.com/](https://deeptime.com/)

---

## **Contact**

For questions or issues, contact Hunter La Force at hunterlaforce61@gmail.com.


